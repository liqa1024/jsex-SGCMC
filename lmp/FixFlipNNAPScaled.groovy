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

/**
 * Add:
 * <pre> {@code
 * fix ID group-ID jse lmp.FixFlipNNAPScaled path/to/jnnpot elem1,elem2,... N ratio seed temp type1:type2,type2:type1,... v_scale
 * } </pre>
 * in lammps in file to use
 */
class FixFlipNNAPScaled extends FixFlipNNAP {
    protected FixFlipNNAPScaled(long aFixPtr, String... aArgs) {
        super(aFixPtr, aArgs)
        if (aArgs.size() > 11) {
            def scaleStr = aArgs[11]
            if (scaleStr.startsWith('v_')) {
                scale = Double.NaN
                scaleStr = scaleStr.substring(2)
                scaleIdx = findVariable(scaleStr)
                if (scaleIdx < 0) throw new IllegalArgumentException("Variable $scaleStr not found")
            } else {
                scale = scaleStr as double
                scaleIdx = -1
            }
        } else {
            scale = 1.0d
            scaleIdx = -1
        }
    }
    @CompileStatic
    @Override double modifyEnergyDiff(double de, int typeOld, int typeNew) {
        return de * (scaleIdx<0 ? scale : computeVariable(scaleIdx))
    }
}
