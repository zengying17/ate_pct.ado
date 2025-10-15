capture program drop ate_pct
program define ate_pct, eclass
    version 14.0
	syntax varlist(min=1), [within truew corr(real 1) control(varlist)]

    // Parse the varlist - these are the group indicator variables (d_i^g)
	
    local groupvars `varlist'
    local G: word count `groupvars'
    // Create a tempvar to indicate sample
    tempvar sample
    qui gen `sample' = e(sample)
    
	tempname beta vcov N r
    matrix `beta' = e(b)
    matrix `vcov' = e(V)
    scalar `N' = e(N)
	scalar `r'=`corr'
	
    // Prepare matrices to hold our estimates
    tempname tau w Sigma_tau Sigma_w
    
    // Extract tau and w from regression results
    //quietly {
        // Retrieve tau values from regression results
        // Assuming tau coefficients are stored in e(b)
    matrix `tau' = J(`G', 1, .)
    matrix `Sigma_tau' = J(`G', `G', 0)
    
    // Extract coefficient estimates and variance-covariance matrix for group dummies
    forvalues g = 1/`G' {
        local dvar: word `g' of `groupvars'
        local col = colnumb(`beta', "`dvar'")
		
        if (`col' > 0) {
            matrix `tau'[`g', 1] = `beta'[1, `col']
            forvalues h = 1/`G' {
                local dvar2: word `h' of `groupvars'
                local col2 = colnumb(`beta', "`dvar2'")
                if (`col2' > 0) {
                    matrix `Sigma_tau'[`g', `h'] = `vcov'[`col', `col2']
                }
            }
        }
        else {
            di as error "Group variable `dvar' not found in regression!"
            exit 198
        }
    }
    

    // Calculate group sizes and weights
    tempvar s_i
    qui gen `s_i' = 0
    foreach dvar of varlist `groupvars' {
    qui    replace `s_i' = `s_i' + `dvar'
    }
    
    tempname N_g N_s p_s w
    matrix `N_g' = J(`G', 1, .)
    
    forvalues g = 1/`G' {
        local dvar: word `g' of `groupvars'
        quietly count if `dvar' == 1 & `sample'
        matrix `N_g'[`g', 1] = r(N)
    }
    
    quietly count if `s_i' == 1 & `sample'
    scalar `N_s' = r(N)
    scalar `p_s' = `N_s'/`N'
    
    // Calculate weights
    matrix `w' = `N_g' * (1/`N_s')
    
    // Calculate weight variance-covariance matrix
    tempname diag_w ww
    matrix `diag_w' = diag(`w')
    matrix `ww' = `w' * `w''
	if ("`truew'" == "") {
    	matrix `Sigma_w' = (1/`N') * (1/`p_s') * (`diag_w' - `ww')
	}
	else {
		matrix `Sigma_w' = J(`G',`G',0)
	}

 
 //  mata: calc_basic()
    
    // If within option is enabled, calculate additional metrics
    if "`within'" != "" {
	local singlecontrol=0
    if "`control'"=="" {
        tempvar uniqctrl
        gen `uniqctrl'=1-`s_i'
		local singlecontrol=1
		
        forvalues j=1/`G' {
            local control "`control' `uniqctrl'"
        }
		}
        else {
            local Gctrl: word count `control'
            if `Gctrl'==1 {
			local singlecontrol=1
				local uniqctrl "`control'"
				forvalues j=2/`G' {
				local control "`control' `uniqctrl'"
				}
			}
            else if `Gctrl'!=`G' {
            di as error "Number of control dummies should either be one or the same as treatment dummies"
            exit 198
            }
        } 		
	tempname N_0
    matrix `N_0' = J(`G', 1, .)    
    forvalues g = 1/`G' {
        local cvar: word `g' of `control'
        quietly count if `cvar' == 1 & `sample'
        matrix `N_0'[`g', 1] = r(N)
    }
	
	 tempname s s1 s0
 quietly {
            // Get residuals from the regression
            tempvar resid residsq
            predict `resid' if `sample', residuals
            gen double `residsq'= `resid'^2
            // Calculate s_g for each group
            matrix `s1' = J(`G', 1, .)
            matrix `s0' = J(`G', 1, .)
            local i = 1
            foreach var of local groupvars {
                summarize `residsq' if `var' == 1 & `sample', meanonly
                matrix `s1'[`i', 1] = sqrt(r(mean))
                local i = `i' + 1
            }
               local i = 1
            foreach var of local control {
                summarize `residsq' if `var' == 1 & `sample', meanonly
                matrix `s0'[`i', 1] = sqrt(r(mean))
                local i = `i' + 1
            }
	        matrix `s'=(`s1' \ `s0')
}
 
mata: 	calc_within("`tau'","`w'","`Sigma_tau'","`Sigma_w'","`s'", "`N_g'", "`N_0'", "`r'", `singlecontrol')

tempname b V
matrix `b'=r(b)
matrix `V'=r(V)

matname `b' taubar rho_a rho_b rho_c rho_b_plus rho_c_plus,explicit c(.)
matname `V' taubar rho_a rho_b rho_c rho_b_plus rho_c_plus,explicit


ereturn post `b' `V'
ereturn local  cmd  "ate_pct"
ereturn scalar p_s = `p_s'

matrix delta=r(delta)
matrix Sigma_delta=r(Sigma_delta)
ereturn matrix delta=delta
ereturn matrix Sigma_delta=Sigma_delta

ereturn display
	}
	else {

		mata:  	calc_basic("`tau'","`w'","`Sigma_tau'","`Sigma_w'")
		tempname b V
matrix `b'=r(b)
matrix `V'=r(V)

matname `b' taubar rho_a rho_b rho_c,explicit c(.)
matname `V' taubar rho_a rho_b rho_c,explicit

ereturn post `b' `V'
ereturn local  cmd  "ate_pct"
ereturn scalar p_s = `p_s'


matrix delta=r(delta)
matrix Sigma_delta=r(Sigma_delta)
ereturn matrix delta=delta
ereturn matrix Sigma_delta=Sigma_delta

ereturn display
	}
end
	
 mata:
 void calc_basic(string scalar tau_name, string scalar w_name, string scalar Sigma_tau_name, string scalar Sigma_w_name)
 {
	tau       = st_matrix(tau_name)
    w         = st_matrix(w_name)
    Sigma_w   = st_matrix(Sigma_w_name)
    Sigma_tau = st_matrix(Sigma_tau_name)
    G = rows(tau)	
    
    eta0      = ln(w) :+ tau
	Sigma_eta = diag(1:/w) * Sigma_w * diag(1:/w) + Sigma_tau
	
	taubar    = w'*tau
	Var_taubar= w'*Sigma_tau*w + tau'*Sigma_w*tau
	rho_a     = exp(taubar)-1
	Var_rho_a = (exp(taubar))^2 * Var_taubar
	
	rho_b     = sum(w :* exp(tau)) - 1
	Var_b     = (exp(eta0))' * Sigma_eta * exp(eta0)
	rho_c=sum(w :* exp(tau :- 0.5 :* diagonal(Sigma_tau))) - 1
	

	b  = (taubar,rho_a,rho_b, rho_c)
	se = (sqrt(Var_taubar),sqrt(Var_rho_a),sqrt(Var_b), sqrt(Var_b))
	
	V  = diag(se:^2)

	st_rclear()
	st_matrix("r(b)",b)
	st_matrix("r(V)",V)       
	st_matrix("r(delta)",(tau\w))
	st_matrix("r(Sigma_delta)",blockdiag(Sigma_tau,Sigma_w))   
}

void calc_within(string scalar tau_name, string scalar w_name, string scalar Sigma_tau_name, string scalar Sigma_w_name, string scalar s_name, string scalar N_g_name, string scalar N_0_name, string scalar r_name, real scalar singlecontrol)
{	
   // The first part repeat calc_basic
	tau       = st_matrix(tau_name)
	w         = st_matrix(w_name)
	Sigma_w   = st_matrix(Sigma_w_name)
	Sigma_tau = st_matrix(Sigma_tau_name)
	taubar    = w'*tau
	Var_taubar= w'*Sigma_tau*w + tau'*Sigma_w*tau
	rho_a     = exp(taubar)-1
	Var_rho_a = (exp(taubar))^2 * Var_taubar 
	
	eta0      = ln(w) :+ tau
	Sigma_eta = diag(1:/w) * Sigma_w * diag(1:/w) + Sigma_tau

	rho_b     = sum(w :* exp(tau)) - 1
	Var_b     = (exp(eta0))' * Sigma_eta * exp(eta0)
	rho_c     = sum(w :* exp(tau :- 0.5 :* diagonal(Sigma_tau))) - 1
	
// Within-group heterogeneity
	s   = st_matrix(s_name)
	N_g = st_matrix(N_g_name)
	N_0 = st_matrix(N_0_name)
	r   = st_numscalar(r_name)
    G   = rows(tau)

	s_g = s[1..G]
	s_0 = s[(G+1)..(2*G)]
	// Calculate Vartau
	Vartau = s_g:^2 + s_0:^2 - 2*r*s_g:*s_0
	// Calculate eta^+ and eta0^+
	eta = ln(w) + tau - 0.5 * (diagonal(Sigma_w) :/ (w:^2)) - 0.5 * diagonal(Sigma_tau)
	eta_plus = eta + 0.5 * Vartau
	eta0_plus = eta0 + 0.5 * Vartau	
	// Calculate rho_b^+ and rho_c^+
	rho_b_plus = sum(w :* exp(tau + 0.5 * Vartau))-1
	rho_c_plus = sum(w :* exp(tau + 0.5 * Vartau - 0.5 * diagonal(Sigma_tau)))-1
	
	// Calculate Sigma_s matrix
	Sigma_s1 = J(G, G, 0)
	Sigma_s0 = J(G, G, 0)
	// Diagonal elements
	for (i=1; i<=G; i++) {
		Sigma_s1[i,i] = 2/(N_g[i]) * s_g[i]^4
        Sigma_s0[i,i] = 2/(N_0[i]) * s_0[i]^4
	}
	
	Sigma_s = blockdiag(Sigma_s1,Sigma_s0)
	Sigma_delta=blockdiag(Sigma_tau,blockdiag(Sigma_w,Sigma_s))

	// Calculate nabla_s0_varsigma
	nabla_sg_varsigma = J(G, 1, 0)
    nabla_s0_varsigma = J(G, 1, 0)
	for (i=1; i<=G; i++) {
		nabla_s0_varsigma[i] = 1 - (r*s_g[i])/s_0[i]
		nabla_sg_varsigma[i] = 1 - (r*s_0[i])/s_g[i]
	}
	
	// Calculate nabla_delta_rho_plus
	
	if (singlecontrol==1){
	nabla_delta_rho_plus = (exp(eta_plus) \ exp(eta_plus):/w \ 0.5 * nabla_sg_varsigma:* exp(eta_plus) \ sum(0.5*nabla_s0_varsigma:*exp(eta_plus)))
	Sigma_delta=Sigma_delta[1..(3*G+1),1..(3*G+1)]
	}
	else{
	nabla_delta_rho_plus = (exp(eta_plus) \ exp(eta_plus):/w \ 0.5 * nabla_sg_varsigma:* exp(eta_plus) \ 0.5*nabla_s0_varsigma:*exp(eta_plus))
	}

	// Calculate Var_b_plus
	Var_b_plus = nabla_delta_rho_plus' * Sigma_delta * nabla_delta_rho_plus
	// Create coefficient and variance matrices
	b = (taubar,rho_a, rho_b, rho_c, rho_b_plus, rho_c_plus)
	se = (sqrt(Var_taubar),sqrt(Var_rho_a),sqrt(Var_b), sqrt(Var_b), sqrt(Var_b_plus), sqrt(Var_b_plus))
	V = diag(se:^2)
	
st_rclear()
st_matrix("r(b)",b)
st_matrix("r(V)",V)    
st_matrix("r(delta)",( tau \ w \ s:^2))
st_matrix("r(Sigma_delta)",Sigma_delta)  
}
end
