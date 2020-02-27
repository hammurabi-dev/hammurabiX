# hammurabi X online documentation

- [installation instruction](./index/install.md)

- [version info](./index/version.md)

- [design info](./index/design.md)

- [performance & precision report](./index/perf.md)

- [author list](./index/author.md)

### about contribution:

We welcome contributions to the source code and documentation from all, with simple rules:

- Please **fork** from the **fix** branch first and create you **pull request** to the **fix** branch, except for some simple modifications.

- Notice that the C++ source code follows certain **indentation** style, and there exist a simple script for auto-correct that, for example, after adding your modifications to your local git, try
```
$ ./contrib/utilities/indent-all
```
this script will try to download some extra files (which is excluded from the git by .gitignore) and apply the indentation automatically, what you need to do is to double check the modifications by
```
$ git status
```
(We acknowledge indentation scripts copied from the ``deal.II`` library.)

- Please let us know if you would like to join us as a contributer or co-author in the [author list](./index/author.md).

- In terms of scientific suggestions or if you want to enrich the physical models or numerical algorithms, you are very welcomed to contact the [authors](./index/author.md) who are currently taking care of the development.
