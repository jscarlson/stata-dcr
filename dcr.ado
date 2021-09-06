*! dcr v5.0.0 JSCarlson 28july2021

***********************
**                   **
**       dcr         **
**                   **
***********************
* Author: Jacob Carlson
* Email: jacob.carlson.research@gmail.com

capture program drop dcr

program dcr, eclass

version 15
syntax anything(id="command line" name=command_line) [if] [in] [fweight aweight pweight iweight], DM1(varname) DM2(varname) DOFUndo(string) [XLSXPATH(string) RSE CRSE(string) SYNTAX(string) PREFIX(string) PSD CHARdm DOFCorr DYADID(string) *]

* -----------------------------------------------
* Check for sortrows
* -----------------------------------------------

capture findfile sortrows.ado

if "`r(fn)'" == "" {
	di as txt "User-written package sortrows needs to be installed first;"
	di as txt "use -ssc install sortrows- to do that."
	exit 498
}

* -----------------------------------------------
* Create options
* -----------------------------------------------

if ("`syntax'" == "") local syntax "vce"
if ("`dofundo'" == "") local dofundo "none"
local remove crse(`crse') syntax(`syntax') xlsxpath(`xlsxpath') rse psd chardm dofcorr dofundo(`dofundo') prefix(`prefix') dyadid(`dyadid')
local options : list options - remove
local oricmd "`: word 1 of `command_line''"

* -----------------------------------------------
* Make varlist
* -----------------------------------------------

foreach avar of varlist _all {
	local allvars "`allvars'" " " "`avar'"
}

local cmdterms : list command_line - allvars
local varlist : list command_line - cmdterms

* -----------------------------------------------
* Create (undirected) dyad ID
* -----------------------------------------------

if ("`dyadid'" != "") {
	tempvar did
	qui gen `did' = `dyadid'
}
else {
	tempvar sdm1 sdm2 did
	sortrows `dm1' `dm2', gen(`sdm1' `sdm2')
	qui egen `did' = group(`sdm1' `sdm2')
}

* -----------------------------------------------
* Mark sample
* -----------------------------------------------

marksample touse
markout `touse' `varlist' `dm1' `dm2', strok

if substr("`oricmd'", 1, 8) == "heckprob" {
	di "Heckman probit being used; marksample turned off..."
	replace `touse' = 1
}

* -----------------------------------------------
* Weight
* -----------------------------------------------

if ("`weight'" != "") local weight [`weight' `exp']

* -----------------------------------------------
* Estimator black-list 
* -----------------------------------------------

error "`oricmd'"=="clogit"

* -----------------------------------------------
* Get list of unique dyad members
* -----------------------------------------------

qui levelsof(`dm1') if `touse', local(unique_dyadmem1)
qui levelsof(`dm2') if `touse', local(unique_dyadmem2)
local unique_dyadmem: list unique_dyadmem1 | unique_dyadmem2

* -----------------------------------------------
* Get number of unique dyad members
* -----------------------------------------------

local N_ud : list sizeof local(unique_dyadmem)
scalar N_ud = `N_ud'

* -----------------------------------------------
* Create contains dyad member clustering variables
* -----------------------------------------------

foreach v in `unique_dyadmem' {
	tempvar contains_`v'
	qui gen `contains_`v'' = .
	if "`chardm'" == "" {
		qui replace `contains_`v'' = -99 if `v'==`dm1' | `v'==`dm2'
		qui replace `contains_`v'' = _n if `v'!=`dm1' & `v'!=`dm2'
	}
	else {
		qui replace `contains_`v'' = -99 if "`v'"==`dm1' | "`v'"==`dm2'
		qui replace `contains_`v'' = _n if "`v'"!=`dm1' & "`v'"!=`dm2'
	}
}

* -----------------------------------------------
* Calculate repeated dyad clustered vcov
* -----------------------------------------------

if "`syntax'" == "vce" {
	qui `prefix' `command_line' if `touse' `weight', `options' vce(cluster `did')
}
else if "`syntax'" == "clusterOnly" {
	qui `prefix' `command_line' if `touse' `weight', `options' cluster(`did')
}
else if "`syntax'" == "robustCluster" {
	qui `prefix' `command_line' if `touse' `weight', `options' cluster(`did')
}

_ms_omit_info e(b)
scalar N_omit = r(k_omit)
matrix temp_V = e(V)
scalar N_par = colsof(temp_V) - N_omit
scalar N_did = e(N_clust)
if N_did == . {
	qui unique `did' if `touse'
	scalar N_did = r(unique)
}
scalar N_obs = e(N)
matrix vrep = e(V)

if "`dofundo'"=="reglike" {
	scalar unadjust = ((N_did - 1)/N_did)*((N_obs - N_par)/(N_obs - 1))
}
else if "`dofundo'"=="none" {
	scalar unadjust = 1
}
else if "`dofundo'"=="asylike" {
	scalar unadjust = ((N_did - 1)/N_did)
}

matrix vcov = -1*unadjust*vrep

* -----------------------------------------------
* Calculate heteroskedasticity robust vcov
* -----------------------------------------------

if "`syntax'" == "vce" {
	qui `prefix' `command_line' if `touse' `weight', `options' vce(robust)
}
else if "`syntax'" == "clusterOnly" {
	qui `prefix' `command_line' if `touse' `weight', `options'
}
else if "`syntax'" == "robustCluster" {
	qui `prefix' `command_line' if `touse' `weight', `options' robust
}

_ms_omit_info e(b)
scalar N_omit = r(k_omit)
scalar N_obs = e(N)
matrix temp_V = e(V)
scalar N_par = colsof(temp_V) - N_omit
matrix vhet = e(V)

if "`dofundo'"=="reglike" {
	scalar unadjust = (N_obs - N_par)/(N_obs)
}
else if "`dofundo'"=="none" {
	scalar unadjust = 1
}
else if "`dofundo'"=="asylike" {  
	scalar unadjust = (N_obs - 1)/(N_obs)
} 

matrix vcov = vcov - ((N_ud - 2)*unadjust*vhet)

* ----------------------------------------------- 
* Calculate contains dyad member clustered vcov's
* -----------------------------------------------

di
_dots 0, title("Multiway decomposition computation in progress!") reps(`N_ud')

foreach v in `unique_dyadmem' {

	if "`syntax'" == "vce" {
		qui `prefix' `command_line' if `touse' `weight', `options' vce(cluster `contains_`v'')
	}
	else if "`syntax'" == "clusterOnly" {
		qui `prefix' `command_line' if `touse' `weight', `options' cluster(`contains_`v'')
	}
	else if "`syntax'" == "robustCluster" {
		qui `prefix' `command_line' if `touse' `weight', `options' cluster(`contains_`v'')
	}

	matrix dspecv = e(V)

	* loading bar
	local i = `i'+1
	_dots `i' 0

	_ms_omit_info e(b)
	scalar N_omit = r(k_omit)
	matrix temp_V = e(V)
	scalar N_par = colsof(temp_V) - N_omit
	scalar N_dyad_`v' = e(N_clust)
	if N_dyad_`v' == . {
    		qui unique `contains_`v'' if `touse'
		scalar N_dyad_`v' = r(unique)
	}
	scalar N_obs = e(N)

	if "`dofundo'"=="reglike" {
		scalar unadjust = ((N_dyad_`v' - 1)/N_dyad_`v')*((N_obs - N_par)/(N_obs - 1))
	}
	else if "`dofundo'"=="none" {
		scalar unadjust = 1
	}
	else if "`dofundo'"=="asylike" {
		scalar unadjust = ((N_dyad_`v' - 1)/N_dyad_`v')
	}

	matrix vcov = (unadjust*dspecv) + vcov

}

* -----------------------------------------------
* Eigendecomposition (code cited from reghdfe)
* -----------------------------------------------

if ("`psd'"!="") {

	tempname Eigenvec lambda lambdaN lambdap lambdapdiag v1N detv1
	
	local v1N = colsof(vcov)
	
	scalar Indnegvar = 0
	
	forvalues l=1(1)`v1N' {
		if vcov[`l',`l']<0 {
			scalar Indnegvar = 1
		}
	}
	
	/*
	matrix `detv1'=det(vcov)
	scalar Inddet = 0
	if `detv1'[1,1]<1e-10 {
		scalar Inddet = 1
	}
	*/
	
	if (Indnegvar==1) {
		matrix symeigen `Eigenvec' `lambda' = vcov
		local lambdaN = colsof(`lambda')
		matrix `lambdap' = `lambda'
		forvalues l=1(1)`lambdaN' {
			if `lambda'[1,`l']>0 {
				matrix `lambdap'[1,`l'] = `lambda'[1,`l']
			}
			else {
				matrix `lambdap'[1,`l'] = 0
			}
		}
		matrix `lambdapdiag' = diag(`lambdap')
		matrix vcov =`Eigenvec'*`lambdapdiag'*(`Eigenvec')'
	}	

}

* -----------------------------------------------
* Apply small sample DOF correction
* -----------------------------------------------

if ("`dofcorr'"!="") {
	matrix vcov = (N_ud/(N_ud - 1)) * vcov
}

* -----------------------------------------------
* Get standard errors for excel output
* -----------------------------------------------

if ("`xlsxpath'"!="") {

	mata: st_matrix("dcr_se", sqrt(diagonal(st_matrix("vcov"))))

	qui putexcel set "`xlsxpath'", replace
	qui putexcel D1="DCR_SE"
	qui putexcel D2=matrix(dcr_se)

	***************************
	* Recreate original model
	***************************

	if ("`crse'"!="") {
		if "`syntax'" == "vce" {
			qui `prefix' `command_line' if `touse' `weight', `options' vce(cluster `crse')
		}
		if "`syntax'" == "clusterOnly" {
			qui `prefix' `command_line' if `touse' `weight', `options' cluster(`crse')
		}
		if "`syntax'" == "robustCluster" {
			qui `prefix' `command_line' if `touse' `weight', `options' cluster(`crse')
		}
	}
	else if ("`crse'"!="") {
		if "`syntax'" == "vce" {
			qui `prefix' `command_line' if `touse' `weight', `options' vce(robust)
		}
		if "`syntax'" == "clusterOnly" {
			qui `prefix' `command_line' if `touse' `weight', `options'
		}
		if "`syntax'" == "robustCluster" {
			qui `prefix' `command_line' if `touse' `weight', `options' robust
		}
	}
	else {
		qui `prefix' `command_line' if `touse' `weight', `options'
	}

	***************************
	* Get info about original clusters
	***************************

	
	if ("`crse'"!="") {
		scalar N_clu_ori = e(N_clust)
		if N_clu_ori == . {
			qui unique `crse' if `touse'
			scalar N_clu_ori = r(unique)
		}
	}
	else {
		scalar N_clu_ori = e(N)
	}

	***************************
	* Get analogous dof correction 
	***************************

	if ("`dofcorr'"!="") {
		matrix ori_V = e(V) * (N_clu_ori/(N_clu_ori - 1))
	} 
	else {
		matrix ori_V = e(V)
	}

	***************************
	* Get parameters
	***************************

	matrix namemat = e(b)'
	local rn: rownames namemat
	local j = 1
	foreach par in `rn' {
		local j = `j' + 1
		qui putexcel A`j'="`par'"
	}

	***************************
	* Get parameter estimates
	***************************

	qui putexcel B1 = "PAR_EST"
	qui putexcel B2 = matrix(e(b)')

	***************************
	* Get original SE
	***************************

	mata: st_matrix("ori_se", sqrt(diagonal(st_matrix("ori_V"))))
	qui putexcel C1 = "ORI_SE"
	qui putexcel C2 = matrix(ori_se)

	***************************
	* Perform SER calculation
	***************************

	mata: st_matrix("ser", st_matrix("dcr_se") :/ st_matrix("ori_se"))
	qui putexcel E1="RATIO_SE"
	qui putexcel E2=matrix(ser)

	***************************
	* P-value (ORI)
	***************************

	mata: st_matrix("ori_pvz", 2*(1 :- normal(abs(st_matrix("e(b)")' :/ st_matrix("ori_se")))))
	qui putexcel F1="ORI_P_Z"
	qui putexcel F2=matrix(ori_pvz)

	if ("`dofcorr'"!="") {
		* P-value from t (ORI)
		local k: list sizeof local(rn)
		scalar k = `k'
		scalar n = e(N)
		mata: st_matrix("ori_pvt", 2*(ttail(st_numscalar("N_clu_ori") - 1, abs(st_matrix("e(b)")' :/ st_matrix("ori_se")))))
		qui putexcel G1="ORI_P_T"
		qui putexcel G2=matrix(ori_pvt)
	} 
	else {
		* P-value from t (ORI)
		local k: list sizeof local(rn)
		scalar k = `k'
		scalar n = e(N)
		mata: st_matrix("ori_pvt", 2*(ttail(st_numscalar("n") - st_numscalar("k"), abs(st_matrix("e(b)")' :/ st_matrix("ori_se")))))
		qui putexcel G1="ORI_P_T"
		qui putexcel G2=matrix(ori_pvt)
	}

	***************************
	* P-value (DCR)
	***************************

	mata: st_matrix("dcr_pv", 2*(1 :- normal(abs(st_matrix("e(b)")' :/ st_matrix("dcr_se")))))
	qui putexcel H1="DCR_P_Z"
	qui putexcel H2=matrix(dcr_pv)

	if ("`dofcorr'"!="") {
		* P-value from t (DCR)
		mata: st_matrix("dcr_pv2", 2*(ttail(st_numscalar("N_ud") - 1, abs(st_matrix("e(b)")' :/ st_matrix("dcr_se")))))
		qui putexcel I1="DCR_P_T"
		qui putexcel I2=matrix(dcr_pv2)

	} 
	else {
		* P-value from t (DCR)
		mata: st_matrix("dcr_pv2", 2*(ttail(st_numscalar("n") - st_numscalar("k"), abs(st_matrix("e(b)")' :/ st_matrix("dcr_se")))))
		qui putexcel I1="DCR_P_T"
		qui putexcel I2=matrix(dcr_pv2)
	}

}

* -----------------------------------------------
* ereturn
* -----------------------------------------------

di
di

matrix oldb = e(b)
scalar numobs = e(N)

ereturn clear
if ("`dofcorr'"!="") { 
	local ssdof = N_ud - 1
	ereturn post oldb vcov, dof(`ssdof')
} 
else {
	ereturn post oldb vcov
}
ereturn scalar N = numobs
ereturn scalar N_unique_dm = N_ud
ereturn display

end
