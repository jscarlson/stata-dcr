program nclogit, eclass
	version 16.1
	syntax varlist [if] [in] [pweight] [, /*
		*/ ROBUST STRata(varname) CLuster(varname) ]
	marksample touse
	markout `touse' `strata' `cluster', strok
	if "`weight'"!="" {
		tempvar w
		quietly generate double `w' `exp' if `touse'
		local iwexp "[iw=`w']"
		capture assert `w' >= 0 if `touse'
		if c(rc) error 402
	}
	if "`cluster'"!="" {
		local clopt "cluster(`cluster')"
	}
	if "`strata'"!="" {
		local stopt "strata(`strata')"
	}

	if "`robust'"!="" | "`cluster'"!="" {
		tempvar s
		tempname D b
		quietly {
			clogit `varlist' `iwexp' if `touse', `stopt'
			matrix `D' = e(V)
			matrix `b' = e(b)
			predict double `s' if e(sample), score
			_robust `s' `iwexp' if e(sample), v(`D') `clopt'
			local n = r(N)
			local n_clust = r(N_clust)
			// return list
			replace `touse' = e(sample)
			ereturn post `b' `D', esample(`touse')
			ereturn scalar N = `n'
			ereturn scalar N_clust = `n_clust'
			ereturn local vcetype "User-Specified"
		}
	}
	else {
		quietly clogit `varlist' `iwexp' if `touse', `stopt'
	}
	
	display
	ereturn display
end
