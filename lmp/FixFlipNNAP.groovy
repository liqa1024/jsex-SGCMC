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
import jse.code.UT
import jsex.nnap.NNAP

/**
 * Add:
 * <pre> {@code
 * fix ID group-ID jse lmp.FixFlipNNAP path/to/jnnpot elem1,elem2,... N ratio seed temp type1:type2,type2:type1,...
 * } </pre>
 * in lammps in file to use
 */
@SuppressWarnings('GrFinalVariableAccess')
class FixFlipNNAP extends FixSgcmcNNAP {
    final Map<Integer, Integer> flipMap
    
    protected FixFlipNNAP(long aFixPtr, String... aArgs) {
        super(aFixPtr)
        
        // parse args
        nnap = new NNAP(aArgs[4])
        def nnapUnits = nnap.units()
        if (nnapUnits) {
            def lmpUnits = unitStyle()
            if (lmpUnits && lmpUnits!=nnapUnits) throw new IllegalArgumentException("Invalid units ($lmpUnits) for this model ($nnapUnits)")
        }
        def elems = aArgs[5].split(',')
        if (elems.size() != ntypes) throw new IllegalArgumentException("Invalid elems size: $elems (ntypes: $ntypes)")
        for (type in 1..elems.size()) {
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
        random = UT.Par.splitRandom(commWorld(), aArgs[8] as long)
        def temps = aArgs[9].split(',')
        if (temps.size() == 1) {
            temp1 = temp2 = temps[0] as double
        } else {
            temp1 = temps[0] as double
            temp2 = temps[1] as double
        }
        flipMap = [:]
        if (aArgs.size() > 10) {
            def tokens = aArgs[10].split(',')
            for (token in tokens) {
                def subArgs = token.split(':')
                int target = subArgs[1] as int
                if (target > ntypes) throw new IllegalArgumentException("Invalid target flip: $target (ntypes: $ntypes)")
                flipMap[subArgs[0] as int] = target
            }
        } else {
            flipMap[1] = 2
            flipMap[2] = 1
        }
    }
    
    @CompileStatic
    @Override int newType(int typeOld) {
        def typeNew = flipMap[typeOld]
        if (typeNew == null) return -1
        return typeNew
    }
    @CompileStatic
    @Override double modifyEnergyDiff(double de, int typeOld, int typeNew) {
        return de
    }
}
