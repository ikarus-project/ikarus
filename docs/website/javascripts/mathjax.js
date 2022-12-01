// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
//
// SPDX-License-Identifier: CC-BY-SA-4.0

window.MathJax = {
  loader: {load: ['[tex]/bm']},
  tex: {
    packages: {'[+]': {'bm'}},
    inlineMath: [['$', '$'], ['\\(', '\\)']],
    displayMath: [['\\[', '\\]']],
    processEscapes: true,
    processEnvironments: true
  },
  options: {ignoreHtmlClass: '.*|', processHtmlClass: 'arithmatex'}
};

document$.subscribe(() => {MathJax.typesetPromise()})