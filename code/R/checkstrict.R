# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  names <- names(found)[found]
  
  if ((length(names) > 0)) {
    sum.nfncs <- 0
    for (i in 1:length(names)) {
      if(!is.function(eval(parse(text=names[i])))) {sum.nfncs <- sum.nfncs + 1}
    }
    if (sum.nfncs > 0) {
      warning("global variables used: ", paste(names(found)[found], collapse=', '))
      return(invisible(FALSE))
    }
  }
  
  !any(found)
}

checkStrict(updateBeta)
checkStrict(updateXi)
checkStrict(updateXiBeta)
checkStrict(updateA)
