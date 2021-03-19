(TeX-add-style-hook "doc"
 (function
  (lambda ()
    (LaTeX-add-bibliographies
     "references")
    (LaTeX-add-labels
     "sec:intro"
     "sec:quick"
     "sec:dataio"
     "sec:algos"
     "sec:exp"
     "sec:plot"
     "sec:data"
     "sec:artdata"
     "sec:realdata")
    (TeX-run-style-hooks
     "latex2e"
     "art10"
     "article"
     "10pt"))))

