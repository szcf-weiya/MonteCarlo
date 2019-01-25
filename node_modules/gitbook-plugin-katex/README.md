Math typesetting using KaTex
==============

Use it for your book, by adding to your book.json:

```
{
    "plugins": ["katex@git+https://github.com/gaoxiaosong/plugin-katex.git"]
}
```

then run `gitbook install`.

## Usage

```
Inline math: $$\int_{-\infty}^\infty g(x) dx$$, $\int_{-\infty}^\infty g(x) dx$


Block math:

$$
\int_{-\infty}^\infty g(x) dx
$$

Or using the templating syntax:

{% block math %}\int_{-\infty}^\infty g(x) dx{% endblock %}
```


### Comparison with [MathJax](https://github.com/GitbookIO/plugin-mathjax)

- Faster

