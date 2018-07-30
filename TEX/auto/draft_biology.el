(TeX-add-style-hook
 "draft_biology"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "english" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "latin9") ("geometry" "margin=1.2in")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "fontenc"
    "inputenc"
    "verbatim"
    "float"
    "amsthm"
    "amsmath"
    "amssymb"
    "graphicx"
    "color"
    "url"
    "caption"
    "subcaption"
    "mathtools"
    "geometry"
    "babel")
   (TeX-add-symbols
    '("inner" 1)
    '("TODO" 1)
    "LL"
    "E"
    "I"
    "ep"
    "Z"
    "GCD"
    "XX"
    "SUM"
    "rr"
    "ii"
    "jj"
    "II"
    "kk"
    "RR"
    "mb"
    "mk"
    "mc"
    "Bell"
    "claimname"
    "definitionname"
    "lemmaname"
    "remarkname"
    "theoremname"
    "corollaryname"
    "propositionname"
    "reals"
    "RL"
    "tamir"
    "CL"
    "RN"
    "RNN"
    "RPP"
    "CNN"
    "hx"
    "one"
    "SNR"
    "be"
    "ee")
   (LaTeX-add-labels
    "fig:micro_example"
    "sec:results"
    "fig:Einst_example"
    "fig:error_per_micro"
    "sec:methods"
    "eq:model"
    "eq:spacing"
    "eq:a2"
    "eq:volume_expansion"
    "eq:projection_model"
    "eq:micrograph_model"
    "eq:Kth_autocorrelation"
    "eq:ac_micrographs"
    "sec:steering"
    "eq:PSWF_defn_eq")
   (LaTeX-add-environments
    "thm"
    "defn"
    "claim"
    "lem"
    "rem"
    "corollary"
    "proposition")
   (LaTeX-add-bibliographies
    "ref"))
 :latex)

