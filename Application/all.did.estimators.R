all.did.estimators <- function(y1, y0, d, x){
  #------------------------------------------------------------------------------------------------------
  #Two-way Fixed Effets
  #twfe.did <- did.twfe.panel(y1, y0, d, x, boot = boot, nboot = nboot)
  twfe.did <- DRDID::twfe_did_panel(y1, y0, d, x, boot = F, nboot = NULL)
  #Compute Regression estimates
  #reg.did <- did.reg.panel.linear(y1, y0, d, x, boot = boot, nboot = nboot)
  reg.did <- DRDID::reg_did_panel(y1, y0, d, x, boot = F, nboot = NULL)
  #Abadie's IPW
  #ipw.did <- did.ipw.panel(y1, y0, d, x, boot = boot, nboot = nboot)
  ipw.did <- DRDID::ipw_did_panel(y1, y0, d, x, boot = F, nboot = NULL)
  #Standardized IPW
  std.ipw.did <- DRDID::std_ipw_did_panel(y1, y0, d, x, boot = F, nboot = NULL)
  # "Traditional" Doubly-Robust: based on OLS regression and MLE pscore
  #dr.trad.did <- dr.did.panel.tr.linear(y1, y0, d, x, x, boot = boot, nboot = nboot)
  dr.trad.did <- DRDID::drdid_panel(y1, y0, d, x, boot = F, nboot = NULL)
  # "Bias-reduced" Doubly-Robust: based on prperly weighted OLS and IPS pscore
  #dr.br.did <- dr.did.panel.br.linear(y1, y0, d, x, boot = boot, nboot = nboot)
  dr.imp.did <- DRDID::drdid_imp_panel(y1, y0, d, x, boot = F, nboot = NULL)
  #------------------------------------------------------------------------------------------------------
  # Put all outputs in a list
  out <- list(twfe = twfe.did,
              reg = reg.did,
              ipw = ipw.did,
              std.ipw = std.ipw.did,
              dr.trad = dr.trad.did,
              dr.imp = dr.imp.did
  )
  #------------------------------------------------------------------------------------------------------
  return(out)
  #------------------------------------------------------------------------------------------------------
}
