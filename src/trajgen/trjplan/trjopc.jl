# Computational Trajectory Generation (CTrjGen)
# Copyright (C) 2025 Qingtan Zeng
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

abstract type AbstTrjOpc end

mutable struct AutoTrjOpc <: AbstTrjOpc

    # 2.0 Constraints
    # 2.1 Dynami system
    dynmdl::DynMdl
    # 2.2, 2.3 states' Ampl-Box and control Ampl&Slop constraints without extra need
    dyncstr::DynCstr
    #2.3 boundaries
    initbrdy::DynBdry
    termbrdy::DynBdry

    #2.4 extra states constraits
    # 2.4.1 collision-free states constraints


end