# Time-Space-Frac-Comp
The code uploaded is a modified version of the lk-wave toolbox, which additionally makes use of the chebyfun toolbox, both avaliable in MATLAB and have their repositiries linked below.

Included modified files from their base packages are;
BirkSongpts.m
  Modified from JACPTS, Gauss-Jacobi quadrature nodes and weights.
  This code has been modified to the best of my knowledge within the copywrite notice given by https://github.com/chebfun/chebfun/blob/master/LICENSE.txt
  The original chebfun repository can be found at https://github.com/chebfun/chebfun
    Copyright (c) 2017, The Chancellor, Masters and Scholars of the University 
    of Oxford, and the Chebfun Developers. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of the University of Oxford nor the names of its 
        contributors may be used to endorse or promote products derived from 
        this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR 
  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

`checkstability.m` `kspaceFirstOrder1DTF.m` `kspaceFirstOrder2DTF.m` `kspaceFirstOrder3DTF.m` `kspaceFirstOrder3DTFv2.m` `kspaceFirstOrder1Dalpha.m` `kspaceFirstOrder2Dalpha.m` `kspaceFirstOrder3Dalpha.m` `private/kspaceFirstOrder_inputChecking.m` `/kspaceFirstOrder_createAbsorptionVariables.m` `/kspaceFirstOrder_expandGridMatricies.m` Modified from files of alike name, and `kspaceFirstOrder1D.m` , `kspaceFirstOrder2D.m` AND `kspaceFirstOrder3D.m` respectively from the kwave toolbox. These code modifications have been completed to the best of the authors knowledge within the GNU Lesser General Public License. Under the initial copywrites Copyright (C) 2009-2019 Bradley Treeby and Ben Cox For further information see http://www.k-wave.org , https://github.com/ucl-bug/k-wave


Each of `kspaceFirstOrder1DTF.m` , `kspaceFirstOrder2DTF.m` , `kspaceFirstOrder3DTF.m` , `kspaceFirstOrder3DTFv2.m` `kspaceFirstOrder1Dalpha.m` , `kspaceFirstOrder2Dalpha.m` and `kspaceFirstOrder3Dalpha.m` are run identically to `kspaceFirstOrder1D.m` etc but `...DTF.m` files require the definition of a scalar `medium.alpha_L` the nuber of Birk-Song quadrature points and weights to be considered. All of these files can run with `medium.alpha_power` and `medium.alpha_coeff` given as scalar or grid-sized variable.
