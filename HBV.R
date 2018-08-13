
# hbvStates = list('sowat' = numeric(), # Soil water storage
#                  'sdep'  = numeric(), # Snow store
#                  'ldep'  = numeric(), # Depth of liquid in snow store
#                  'stw1'  = numeric(), # Shallow soil storage
#                  'stw2'  = numeric()  # Deep soil storage
#                  )

# hbvFluxes = list('qRouting' = numeric(),
#                  'qSim'     = numeric(),
#                  'actualET' = numeric()
#                  )

hbv.sim <- function(data,
                    hl1, # Storage of shallow layer (mm)
                    ck0, # Rate constants
                    ck1,
                    ck2, 
                    perc, # Percolation rate (mm/d)
                    lp, # Evap constant (unitless)
                    fcap, # Surface soil moisture (mm)
                    beta, # Soil storage exponent (unitless)
                    maxbas, # Number of days for hydrograph routing
                    ttlim, # Snow threshold (ºC)
                    degd, # Degree day factor (mm/(ºC-d))
                    degw,
                    return_state = FALSE)
{

  # Check data has been entered
  stopifnot(c('P', 'E', 'T') %in% colnames(data))
  
  inAttr <- attributes(data[,1])
  data <- as.ts(data)
  
  P <- data[,'P']
  Q <- data[,'Q']
  E <- data[,'E']
  airT <- data[,'T']
  
  # Skip missing values
  bad <- is.na(P) | is.na(E) | is.na(airT)
  P[bad] <- 0
  E[bad] <- 0
  airT[bad] <- 0
  
  # Set up vectors
  sowat = rep(0, nrow(data)) # Soil water storage
  sdep  = rep(0, nrow(data)) # Snow store
  ldep  = rep(0, nrow(data)) # Depth of liquid in snow store
  stw1  = rep(0, nrow(data)) # Shallow soil storage
  stw2  = rep(0, nrow(data)) # Deep soil storage
  
  # Fluxes
  qRouting = rep(0, 2)
  qSim     = rep(0, nrow(data))
  actualET = rep(0, nrow(data))
  
  # Model
  for(t in 2:nrow(data)) {
    # ------------ #
    # Snow routine # 
    # ------------ #
    # Reset local pars
    smelt = 0
    effectiveP = 0
    
    # Grab data
    avgTemp <- airT[t]
    precip <- P[t]
    
    # Set starting state equal to previous day
    sdep[t] = sdep[t-1]
    # Snow if below threshold
    if(avgTemp < ttlim) {
      sdep[t] = sdep[t] + precip # Add all P to snow store
    } else {
      effectiveP = precip
    }
    # Snow melt if temperature greater than threshold
    if(avgTemp > degw) {
      # Check there is snow to melt
      if(sdep[t] > 0) {
        smelt = (avgTemp - degw)*degd
        # Check if there is more snow melting than available
        if(smelt > sdep[t]) {
          effectiveP = effectiveP + sdep[t] # Add full snow depth to effective P
          sdep[t] = 0 # All snow melted
        } else {
          # Melt portion of snow
          effectiveP = effectiveP + smelt # effective P + snow melt
          sdep[t] = sdep[t] - smelt # Remove melted snow from store
        }
      }
    }
    # ------------ #
    # Soil routine #
    # ------------ #
    # Local pars
    hsw = 0
    aET = 0
    runoffDepth = 0
    # Grab data
    pET = E[t]
    
    # Starting point equal to yesterday
    sowat[t] = sowat[t-1]
    
    # If soil storage is full, runoff = P + excess
    if(sowat[t] >= fcap) {
      runoffDepth = effectiveP + (sowat[t] - fcap)
      sowat[t] = fcap
    } else {
      # Portion of effectiveP that is stored
      hsw = effectiveP * (1 - (sowat[t]/fcap)^beta)
      sowat[t] = sowat[t] + hsw
      runoffDepth = effectiveP - hsw
      # If soil storage will be exceeded
      if(sowat[t] > fcap) {
        runoffDepth = runoffDepth - (sowat[t] - fcap)
        sowat[t] = fcap # Back to capacity
      }
    }
    
    aET = pET * min(sowat[t-1]/(fcap * lp), 1) # actual ET after adjustment for SM
    if(aET < 0) aET = 0
    
    # If there is soil moisture for aET
    if(sowat[t] > aET) {
      actualET[t] = aET
      sowat[t] = sowat[t] - aET
    } else {
      actualET[t] = sowat[t]
      sowat[t] = 0 # all evaporated
    }
    
    stw1[t] = stw1[t-1] + runoffDepth
    
    # --------- #
    # Discharge #
    # --------- #
    # Local pars
    Q0 = 0
    Q1 = 0
    Q2 = 0
    Qall = 0
    
    # If upper reservoir is above threshold for near surface flow
    if(stw1[t] > hl1) {
      # Calculate and remove from reservoir
      Q0 = (stw1[t] - hl1) * ck0
      stw1[t] = stw1[t] - Q0
    } else {
      Q0 = 0
    }
    
    # If there is water left in the upper reservoir
    if(stw1[t] > 0) {
      # Calculate interflow and remove
      Q1 = stw1[t] * ck1
      stw1[t] = stw1[t] - Q1
    } else {
      Q1 = 0
    }
    
    # If there is enough water to completely supply percolation
    if(stw1[t] > perc) {
      # Move water from upper to lower
      stw1[t] = stw1[t] - perc
      stw2[t] = stw2[t] + perc
    } else {
      # Move what we can
      stw2[t] = stw2[t] + stw1[t]
      stw1[t] = 0
    }
    
    # If there is water in lower reservoir
    if(stw2[t] > 0) {
      # Calculate baseflow and remove it
      Q2 = stw2[t] * ck2
      stw2[t] = stw2[t] - Q2
    } else {
      Q2 = 0 
    }
    
    Qall = Q0 + Q1 + Q2 # All discharge
    
    # ------- #
    # Routing #
    # ------- #
    # Local pars
    m2   = (maxbas / 2) - 1
    wsum = 0
    wei  = rep(0, maxbas)
    
    # wei | g(t,MAXBAS) | transformation function consisting os a triangular weighting function and one free parameter
    # Qrouting | NA | This is the flow from the single Qall spread out over time according to the transformation function
    # Qsim | NA | The final flow output by the model
    
    # Calculate the values of the transformation function according to maxbas
    for(i in 1:maxbas) {
      if(i <= m2) {
        wei[i] = i+1
      } else {
        wei[i] = (maxbas - (i+1)) + 1
      }
      wsum = wsum + wei[i]
    }
    
    # Now, spread the flow Qall out over Qind according to the transformation function
    for(i in 1:maxbas) {
      wei[i] = wei[i]/wsum
      # Qind is constantly added to by the transformed Qall.  In other words, when Qall is transformed (spread out over time)
      # it is then added to whatever currently exists in Qind for those time steps.  In other words, a previous transformation of
      # Qall for the previous time step placed flows in Qind in times that overlapped with the currently transformed flow times.
      qRouting[i] = qRouting[i] + Qall * wei[i]
    }
    
    qSim[t] = qRouting[1]
    
    # -------- #
    # Backflow #
    # -------- #
    klen = 2 * maxbas - 1
    
    for(k in 1:klen) {
      qRouting[k] = qRouting[k+1]
    }
    qRouting[klen] = 0 
  } # close loop
  
  U <- qSim
  U[bad] <- NA # Remove missing values
  
  # Attributes
  attributes(U) <- inAttr
  ans <- U
  
  if (return_state == TRUE) {
    sowat[bad] = NA 
    sdep[bad]  = NA 
    ldep[bad]  = NA 
    stw1[bad]  = NA 
    stw2[bad]  = NA 
    
    qRouting[bad] = NA
    actualET[bad] = NA
    
    attributes(sowat) = inAttr 
    attributes(sdep)  = inAttr 
    attributes(ldep)  = inAttr 
    attributes(stw1)  = inAttr 
    attributes(stw2)  = inAttr 
    
    attributes(qRouting) = inAttr
    attributes(actualET) = inAttr
    

    ans <- cbind(U=U,
                 AcutalET=actualET,
                 SoilWater=sowat,
                 SnowDepth=sdep,
                 LiquidDepth=ldep,
                 UpperS=stw1,
                 LowerS=stw2,
                 qRouting=qRouting)
  }
  # Return ans
  ans
}
