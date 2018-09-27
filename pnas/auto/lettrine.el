(TeX-add-style-hook
 "lettrine"
 (lambda ()
   (TeX-run-style-hooks
    "keyval")
   (TeX-add-symbols
    '("LettrineOptionsFor" 2)
    "DefaultOptionsFile"
    "DefaultLoversize"
    "DefaultLraise"
    "DefaultLhang"
    "DefaultFindent"
    "DefaultNindent"
    "DefaultSlope"
    "L"
    "LettrineTextFont"
    "LettrineFontHook"
    "LettrineFont"
    "LettrineFontEPS"
    "LettrineWidth"
    "color"
    "textcolor")
   (LaTeX-add-counters
    "DefaultLines")
   (LaTeX-add-saveboxes
    "L"))
 :latex)

