```{include} ../README.md
:end-before: '## Citation'
```

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation {octicon}`plug;1em;`
:link: installation
:link-type: doc

New to *gcell*? Check out the installation guide.
:::


:::{grid-item-card} API reference {octicon}`book;1em;`
:link: api/index
:link-type: doc

The API reference contains a detailed description of
the gcell API.
:::


:::{grid-item-card} GitHub {octicon}`mark-github;1em;`
:link: https://github.com/GET-Foundation/gcell

Find a bug? Interested in improving gcell? Checkout our GitHub for the latest developments.
:::
::::

**Other resources**

* Follow changes in the {ref}`release notes <release-notes>`.
* Check out our {ref}`contribution guide <contribution-guide>` for development practices.

### News

```{include} news.md
:start-after: '<!-- marker: after prelude -->'
:end-before: '<!-- marker: before old news -->'
```

{ref}`(past news) <News>`

% put references first so all references are resolved

% NO! there is a particular meaning to this sequence

```{toctree}
:hidden: true
:maxdepth: 1

installation
api/index
release-notes/index
community
news
dev/index
contributors
references
```

[contribution guide]: dev/index.md
[Nature (2024)]: https://www.nature.com/articles/s41586-024-08391-z
[github]: https://github.com/GET-Foundation/gcell
