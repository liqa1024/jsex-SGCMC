/**
 * This file is part of example scripts of FFS in jse
 * Copyright 2025 Qing'an Li
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
package lmp

import groovy.transform.CompileStatic
import jse.cache.IntVectorCache
import jse.cache.MatrixCache
import jse.clib.IntCPointer
import jse.clib.NestedIntCPointer
import jse.code.UT
import jse.code.collection.IntList
import jse.code.random.IRandom
import jse.lmp.LmpPlugin
import jse.math.MathEX
import jse.math.matrix.RowMatrix
import jse.math.vector.IntVector
import jse.parallel.MPI
import jsex.nnap.NNAP

/**
 * The (vc-)sgcmc algorithm suitable for jse nnap adopts almost the same writing method as
 * fix sgcmc in lammpa, but does not use the switching window algorithm.
 * <p>
 * Add:
 * <pre> {@code
 * fix ID group-ID jse lmp.FixSgcmcNNAP path/to/jnnpot elem1,elem2,... nevery ratio temp deltamu,... [variance kappa conc,...] [randseed seed] [scale v_scale] [recenter]
 * } </pre>
 * in lammps in file to use.
 * <p>
 * Reference:
 * <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.184203">
 * Scalable parallel Monte Carlo algorithm for atomistic simulations of precipitation in alloys </a>
 */
@SuppressWarnings('GrFinalVariableAccess')
class FixSgcmcNNAP extends LmpPlugin.Fix {
    final double forceBoltz
    final int ntypes
    NNAP nnap
    final int[] lmpType2nnapType
    final double[] cutoff, cutsq
    double cutoffmax, cutmaxsq, cutoff2, cut2sq
    int nevery_
    double flipRatio
    double temp1, temp2
    final double[] deltamu
    boolean variance = false
    double kappa = 0.0d
    double[] conc = null
    Long seed = null
    IRandom random
    double scale = 1.0d
    int scaleIdx = -1
    boolean recenter = false
    
    long flipSucNum = 0, flipNum = 0
    int[] natomsType, natomsTypeDelta
    
    protected FixSgcmcNNAP(long aFixPtr) {
        super(aFixPtr)
        // basic setting
        setDynamicGroupAllow(true)
        setVectorFlag(true)
        setSizeVector(2)
        setGlobalFreq(1)
        setExtvector(false)
        setTimeDepend(true)
        setForceReneighbor(true)
        setNextReneighbor(ntimestep()+1)
        forceBoltz = forceBoltz()
        ntypes = atomNtypes()
        if (ntypes < 2) throw new IllegalArgumentException('Type number of atom MUST >= 2')
        
        natomsType = new int[ntypes+1]
        natomsTypeDelta = new int[ntypes+1]
        lmpType2nnapType = new int[ntypes+1]
        cutoff = new double[ntypes+1]
        cutsq = new double[ntypes+1]
        deltamu = new double[ntypes+1]
    }
    
    protected FixSgcmcNNAP(long aFixPtr, String... aArgs) {
        this(aFixPtr)
        
        // parse args
        nnap = new NNAP(aArgs[4])
        def nnapUnits = nnap.units()
        if (nnapUnits) {
            def lmpUnits = unitStyle()
            if (lmpUnits && lmpUnits!=nnapUnits) throw new IllegalArgumentException("Invalid units ($lmpUnits) for this model ($nnapUnits)")
        }
        def elems = aArgs[5].split(',')
        if (elems.size() != ntypes) throw new IllegalArgumentException("Invalid elems size: $elems (ntypes: $ntypes)")
        for (type in 1..ntypes) {
            def elem = elems[type-1]
            int nnapType = nnap.typeOf(elem)
            if (nnapType <= 0) throw new IllegalArgumentException("Invalid element ($elem) in pair_coeff")
            lmpType2nnapType[type] = nnapType
            cutoff[type] = nnap.model(nnapType).basis().rcut()
            cutsq[type] = cutoff[type]*cutoff[type]
        }
        cutoffmax = cutoff.max(); cutmaxsq = cutoffmax*cutoffmax
        cutoff2 = cutoffmax+cutoffmax; cut2sq = cutoff2*cutoff2
        nevery_ = aArgs[6] as int
        setNevery(nevery_)
        flipRatio = aArgs[7] as double
        def temps = aArgs[8].split(',')
        if (temps.size() == 1) {
            temp1 = temp2 = temps[0] as double
        } else {
            temp1 = temps[0] as double
            temp2 = temps[1] as double
        }
        def deltamu_ = aArgs[9].split(',')
        if (deltamu_.size()+1!=ntypes) throw new IllegalArgumentException("Invalid deltamu size: $deltamu_ (ntypes: $ntypes)")
        for (type in 2..ntypes) {
            deltamu[type] = deltamu_[type-2] as double
        }
        int idx = 10
        while (idx > 0) idx = parseArgs(aArgs, idx)
        random = seed==null ? UT.Par.splitRandom(commWorld()) : UT.Par.splitRandom(commWorld(), seed)
    }
    
    private int parseArgs(String[] aArgs, int idx) {
        if (aArgs.size() <= idx) return -1
        def key = aArgs[idx]
        if (key == 'variance') {
            if (conc != null) throw new IllegalArgumentException("Redundant key: $key")
            variance = true
            kappa = aArgs[idx+1] as double
            def conc_ = aArgs[idx+2].split(',')
            if (conc_.size()+1!=ntypes) throw new IllegalArgumentException("Invalid conc size: $conc_ (ntypes: $ntypes)")
            conc = new double[ntypes+1]
            conc[1] = 1.0d
            for (type in 2..ntypes) {
                double conci = conc_[type-2] as double
                if (conci < 0) throw new IllegalArgumentException("Concentration MUST NOT be negative")
                conc[type] = conci
                conc[1] -= conci
            }
            if (conc[1] < 0) throw new IllegalArgumentException("Total concentration exceeds 1")
            return idx + 3
        } else
        if (key == 'randseed') {
            if (seed != null) throw new IllegalArgumentException("Redundant key: $key")
            seed = aArgs[idx+1] as long
            return idx + 2
        } else
        if (key == 'scale') {
            def scaleStr = aArgs[idx+1]
            if (scaleStr.startsWith('v_')) {
                scale = Double.NaN
                scaleStr = scaleStr.substring(2)
                scaleIdx = findVariable(scaleStr)
                if (scaleIdx < 0) throw new IllegalArgumentException("Variable $scaleStr not found")
            } else {
                scale = scaleStr as double
                scaleIdx = -1
            }
            return idx + 2
        } else if (key == 'recenter') {
            recenter = true
            return idx + 1
        } else {
            throw new IllegalArgumentException("Invalid key: $key")
        }
    }
    
    long lastStep = -1, firstStep = -1
    @CompileStatic
    @Override void init() {
        lastStep = laststep()
        firstStep = firststep()
        // Neighbor list requiring twice cutoff radius
        double skin = neighborSkin()
        double cutrq = cutoff2 + skin
        if (Math.max(forcePairCutforce()+skin, commCutghostuser()) < cutrq) {
            throw new IllegalStateException("FixSgcmcNNAP need 2 x cutoff, use `comm_modify cutoff $cutrq`")
        }
        neighborRequestOccasionalFull(cutrq)
    }
    @Override int setMask() {
        int mask = 0
        mask |= POST_FORCE
        return mask
    }
    @Override double computeVector(int i) {
        if (i == 0) return flipNum
        if (i == 1) return flipSucNum
        return 0.0d
    }
    
    @CompileStatic
    double getBeta() {
        double temp = temp1 + (temp2-temp1)*(ntimestep()-firstStep)/(lastStep-firstStep)
        return 1.0d / (forceBoltz*temp)
    }
    
    final IntList ilist = new IntList(16)
    final IntList nl2 = new IntList(16)
    final IntList nlc = new IntList(16)
    final List<IntList> nlall = []
    @CompileStatic
    @Override void postForce(int aVFlag) {
        long ntimestep_ = ntimestep()
        if (nextReneighbor() != ntimestep_) return
        
        int flipSucNum_ = 0
        int flipNum_ = 0
        def world = commWorld()
        
        def finalDo = {
            flipSucNum += world.allreduceI(flipSucNum_, MPI.Op.SUM)
            flipNum += world.allreduceI(flipNum_, MPI.Op.SUM)
            setNextReneighbor(ntimestep_ + nevery_)
        }
        
        neighborBuildOne()
        final int nlocal = atomNlocal()
        final int nghost = atomNghost()
        final int nmax = nlocal+nghost
        final def numneigh = listNumneigh()
        final def firstneigh = listFirstneigh()
        final int groupbit_ = groupbit()
        final def mask = atomMask()
        
        // Stat the ilist of the current group
        ilist.clear()
        for (int i = 0; i < nlocal; ++i) {
            if ((mask[i]&groupbit_) != 0) ilist << i
        }
        final int nlist = ilist.size()
        double nflips = flipRatio*nlist*nevery_
        // vcsgcmc requires all processes to use the same number of flips, because each flip requires a total synchronization
        int nflipsMax = world.allreduceI(MathEX.Code.ceil2int(nflips), MPI.Op.MAX)
        if (nflipsMax == 0) {
            finalDo()
            return
        }
        
        def mass = atomMass()
        def v = atomV()
        def x = atomX()
        def type = atomType()
        def xMat = MatrixCache.getMatRow(nmax, 3)
        def vMat = MatrixCache.getMatRow(nlocal, 3)
        def typeVec = IntVectorCache.getVec(nmax)
        x.parse2dest(xMat)
        v.parse2dest(vMat)
        type.parse2dest(typeVec)
        finalDo <<= {
            type.fill(typeVec)
            v.fill(vMat)
            IntVectorCache.returnVec(typeVec)
            MatrixCache.returnMat(vMat)
            MatrixCache.returnMat(xMat)
        }
        // Considering that the mask may change, the natomsType is updated uniformly before each step.
        int natoms = -1
        if (variance) {
            natoms = world.allreduceI(nlist, MPI.Op.SUM)
            for (int type_ = 1; type_ <= ntypes; ++type_) {
                natomsType[type_] = 0
            }
            for (int ii = 0; ii < nlist; ++ii) {
                ++natomsType[typeVec[ilist[ii]]]
            }
            world.allreduce(natomsType, 1, ntypes, MPI.Op.SUM)
        }
        
        // Start the flip
        for (int i = 0; i < nflipsMax; ++i) {
            if (variance) {
                for (int type_ = 1; type_ <= ntypes; ++type_) {
                    natomsTypeDelta[type_] = 0
                }
            }
            int choice = -1
            int typeOld = -1
            // The difference is reflected here through probability of non-execution.
            if (random.nextDouble()*nflipsMax < nflips) {
                // Randomly select atom
                choice = ilist[random.nextInt(nlist)]
                typeOld = typeVec[choice]
                // Randomly select target type
                int typeNew = newType(typeOld)
                if (typeNew > 0) {
                    ++flipNum_
                    // Calculate flip energy difference
                    double de = calEnergyDiffFlip(choice, typeNew, numneigh, firstneigh, nlocal, xMat, typeVec)
                    // Modifying this energy difference is required for a specific ensemble
                    de = modifyEnergyDiff(de, typeOld, typeNew)
                    // Perform MC for the current process
                    if (de<=0.0d || random.nextDouble()<MathEX.Fast.exp(-de*beta)) {
                        // Energy conservation scaling here
                        double scale = MathEX.Fast.sqrt((double)mass[typeOld] / (double)mass[typeNew])
//                        double scale = mass[typeOld] / mass[typeNew]
                        vMat.set(choice, 0, scale*vMat.get(choice, 0))
                        vMat.set(choice, 1, scale*vMat.get(choice, 1))
                        vMat.set(choice, 2, scale*vMat.get(choice, 2))
                        ++flipSucNum_
                        if (variance) {
                            --natomsTypeDelta[typeOld]
                            ++natomsTypeDelta[typeNew]
                        }
                    } else {
                        typeVec[choice] = typeOld
                        choice = -1
                        typeOld = -1
                    }
                } else {
                    choice = -1
                    typeOld = -1
                }
            }
            // vcsgcmc needs to perform global flipping again
            if (variance) {
                // Count the total number of modifications
                world.allreduce(natomsTypeDelta, 1, ntypes, MPI.Op.SUM)
                // Only those who perform flipping need to consider the fallback
                if (choice >= 0) {
                    double da = 0.0d
                    for (int type_ = 1; type_ <= ntypes; ++type_) {
                        da += natomsTypeDelta[type_] * natomsTypeDelta[type_]
                        da += 2*natomsTypeDelta[type_] * (natomsType[type_] - conc[type_]*natoms)
                    }
                    da *= (kappa/natoms)
                    // On failure, fallback to the previous flip
                    if (da>0.0d && random.nextDouble()>=MathEX.Fast.exp(-da)) {
                        // Energy conservation scaling here
                        double scale = MathEX.Fast.sqrt((double)mass[typeVec[choice]] / (double)mass[typeOld])
//                        double scale = mass[typeVec[choice]] / mass[typeOld]
                        vMat.set(choice, 0, scale*vMat.get(choice, 0))
                        vMat.set(choice, 1, scale*vMat.get(choice, 1))
                        vMat.set(choice, 2, scale*vMat.get(choice, 2))
                        --flipSucNum_
                        typeVec[choice] = typeOld
                    }
                }
            }
        }
        // Velocity recenter as needed
        if (recenter) {
            double totMass = 0.0d
            double totMomX = 0.0d
            double totMomY = 0.0d
            double totMomZ = 0.0d
            for (int ii = 0; ii < nlist; ++ii) {
                int i = ilist[ii]
                double massi = mass[typeVec[i]]
                totMass += massi
                totMomX += massi * vMat.get(i, 0)
                totMomY += massi * vMat.get(i, 1)
                totMomZ += massi * vMat.get(i, 2)
            }
            totMass = world.allreduceD(totMass, MPI.Op.SUM)
            totMomX = world.allreduceD(totMomX, MPI.Op.SUM)
            totMomY = world.allreduceD(totMomY, MPI.Op.SUM)
            totMomZ = world.allreduceD(totMomZ, MPI.Op.SUM)
            double vcx = totMomX / totMass
            double vcy = totMomY / totMass
            double vcz = totMomZ / totMass
            for (int ii = 0; ii < nlist; ++ii) {
                int i = ilist[ii]
                vMat.set(i, 0, vMat.get(i, 0) - vcx)
                vMat.set(i, 1, vMat.get(i, 1) - vcy)
                vMat.set(i, 2, vMat.get(i, 2) - vcz)
            }
        }
        finalDo()
    }
    
    @CompileStatic
    double modifyEnergyDiff(double de, int typeOld, int typeNew) {
        // sgcmc needs to consider the difference in chemical potential before and after flipping
        de -= deltamu[typeOld]
        de += deltamu[typeNew]
        return de * (scaleIdx<0 ? scale : computeVariable(scaleIdx))
    }
    @CompileStatic
    int newType(int typeOld) {
        if (ntypes == 2) {
            return typeOld==1 ? 2 : 1
        } else {
            int typeNew = random.nextInt(ntypes-1) + 1
            if (typeNew >= typeOld) ++typeNew
            return typeNew
        }
    }
    @CompileStatic
    double calEnergyDiffFlip(int choice, int typeNew, IntCPointer numneigh, NestedIntCPointer firstneigh, int nlocal, RowMatrix xMat, IntVector typeVec) {
        double xc = xMat.get(choice, 0)
        double yc = xMat.get(choice, 1)
        double zc = xMat.get(choice, 2)
        // Get the neighbor list, this list is twice as large so it needs to be traversed manually
        def jlist = firstneigh[choice]
        int jnum = numneigh[choice]
        // First traverse once to obtain a list of twice the cutoff that meets the requirements, as well as a list of neighbors of the central atom.
        nl2.clear()
        nlc.clear()
        nl2 << choice // need to include self
        nlc << choice // need to include self
        for (int jj = 0; jj < jnum; ++jj) {
            int j = jlist[jj]
            j &= LmpPlugin.NEIGHMASK
            double dx = xMat.get(j, 0) - xc
            double dy = xMat.get(j, 1) - yc
            double dz = xMat.get(j, 2) - zc
            double rsq = dx*dx + dy*dy + dz*dz
            if (rsq < cutmaxsq) nlc << j
            if (rsq < cut2sq) nl2 << j
        }
        // Calculate the total energy before flipping, and by the way, cache all neighbor lists (used to speed up energy calculation after flipping)
        double[] oEng = [0.0d]
        final int nnc = nlc.size()
        nlall.forEach {it.clear()}
        while (nlall.size() < nnc) nlall << new IntList(16)
        nnap.calEnergyPart(nlocal, {initDo, finalDo_, nlDo ->
            initDo?.run(0)
            for (int kk = 0; kk < nnc; ++kk) {
                def nlk = nlall[kk]
                int k = nlc[kk]
                double xk = xMat.get(k, 0)
                double yk = xMat.get(k, 1)
                double zk = xMat.get(k, 2)
                int typek = typeVec[k]
                double cutsqk = cutsq[typek]
                double cutsqkNew = k==choice ? cutsq[typeNew] : cutsqk // When the type changes, the cutoff radius will also change.
                nlDo.run(0, k, lmpType2nnapType[typek]) {rmax, dxyzTypeIdxDo ->
                    final int nn = nl2.size()
                    for (int jj = 0; jj < nn; ++jj) {
                        int j = nl2[jj]
                        if (j == k) continue
                        double dx = xMat.get(j, 0) - xk
                        double dy = xMat.get(j, 1) - yk
                        double dz = xMat.get(j, 2) - zk
                        double rsq = dx*dx + dy*dy + dz*dz
                        // When the type changes, the cutoff radius will also change.
                        if (rsq < cutsqkNew) nlk << j
                        if (rsq < cutsqk) dxyzTypeIdxDo.run(dx, dy, dz, lmpType2nnapType[typeVec[j]], j)
                    }
                }
            }
            finalDo_?.run(0)
        }, {threadID, cIdx, eng ->
            oEng[0] += eng
        })
        // Calculate the total energy after flipping
        typeVec[choice] = typeNew
        double[] nEng = [0.0d]
        nnap.calEnergyPart(nlocal, {initDo, finalDo_, nlDo ->
            initDo?.run(0)
            for (int kk = 0; kk < nnc; ++kk) {
                def nlk = nlall[kk]
                int k = nlc[kk]
                double xk = xMat.get(k, 0)
                double yk = xMat.get(k, 1)
                double zk = xMat.get(k, 2)
                int typek = typeVec[k]
                nlDo.run(0, k, lmpType2nnapType[typek]) {rmax, dxyzTypeIdxDo ->
                    final int nn = nlk.size()
                    for (int jj = 0; jj < nn; ++jj) {
                        int j = nlk[jj]
                        double dx = xMat.get(j, 0) - xk
                        double dy = xMat.get(j, 1) - yk
                        double dz = xMat.get(j, 2) - zk
                        dxyzTypeIdxDo.run(dx, dy, dz, lmpType2nnapType[typeVec[j]], j)
                    }
                }
            }
            finalDo_?.run(0)
        }, {threadID, cIdx, eng ->
            nEng[0] += eng
        })
        return nEng[0] - oEng[0]
    }
}
