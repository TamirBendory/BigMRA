(TeX-add-style-hook
 "needspace"
 (lambda ()
   (TeX-add-symbols
    '("needspace" 1)
    "Needspace")
   (LaTeX-add-environments
    '("mdframed" LaTeX-env-args ["argument"] 0)))
 :latex)

