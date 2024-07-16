# Create tables

out.table.lalonde <- function(sample1, benchmark = NULL){
 
  l1.lin <- cbind(sample1$linear$twfe$ATT,
                  sample1$linear$reg$ATT,
                  sample1$linear$ipw$ATT,
                  sample1$linear$std.ipw$ATT,
                  sample1$linear$dr.trad$ATT,
                  sample1$linear$dr.imp$ATT)
  
  l2.lin <- cbind(sample1$linear$twfe$se,
                  sample1$linear$reg$se,
                  sample1$linear$ipw$se,
                  sample1$linear$std.ipw$se,
                  sample1$linear$dr.trad$se,
                  sample1$linear$dr.imp$se)
  
  
  l4.lin <- cbind(sample1$dw$twfe$ATT,
                  sample1$dw$reg$ATT,
                  sample1$dw$ipw$ATT,
                  sample1$dw$std.ipw$ATT,
                  sample1$dw$dr.trad$ATT,
                  sample1$dw$dr.imp$ATT)
  
  l5.lin <- cbind(sample1$dw$twfe$se,
                  sample1$dw$reg$se,
                  sample1$dw$ipw$se,
                  sample1$dw$std.ipw$se,
                  sample1$dw$dr.trad$se,
                  sample1$dw$dr.imp$se)
  
  
  l7.lin <- cbind(sample1$aug.dw$twfe$ATT,
                  sample1$aug.dw$reg$ATT,
                  sample1$aug.dw$ipw$ATT,
                  sample1$aug.dw$std.ipw$ATT,
                  sample1$aug.dw$dr.trad$ATT,
                  sample1$aug.dw$dr.imp$ATT)
  
  l8.lin <- cbind(sample1$aug.dw$twfe$se,
                  sample1$aug.dw$reg$se,
                  sample1$aug.dw$ipw$se,
                  sample1$aug.dw$std.ipw$se,
                  sample1$aug.dw$dr.trad$se,
                  sample1$aug.dw$dr.imp$se)
  
  if(!is.null(benchmark)){
    l3.lin <- l1.lin/benchmark
    l6.lin <- l4.lin/benchmark
    l9.lin <- l7.lin/benchmark
  } else {
    l3.lin <- l1.lin*0
    l6.lin <- l4.lin*0
    l9.lin <- l7.lin*0
  }
  
  
  table <- rbind(l1.lin, l2.lin, l3.lin,
                 l4.lin, l5.lin, l6.lin,
                 l7.lin, l8.lin, l9.lin)
  
  colnames(table) <- c("TWFE", "REG", "IPW","STD-IPW",  "DR-TR", "DR-imp")
  rownames(table) <- c("lin.att", "lin.se", "lin.ev.bias",
                       "dw.att", "dw.se", "dw.ev.bias",
                       "aug.dw.att", "aug.dw.se", "aug.dw.ev.bias")
  return(table)
  
}