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
