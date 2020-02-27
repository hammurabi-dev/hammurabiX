## hammurabi X Pipeline class

The ``Pipeline`` class is designed for assembling the building blocks of the simulation routine.
The ``Pipeline`` class contains unique pointers for all field and grid classes internally, the assembling functions are responsible for binding these pointers correctly and dynamically.
The constructor of ``Pipeline`` class parses the given parameter-set file and stores all parameters in a ``Param`` class object.

- [header file](https://github.com/hammurabi-dev/hammurabiX/tree/master/include/pipeline.h)
- [source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/pipeline/pipeline.cc)

## function list:

- **``Pipeline::Pipeline``**
```
# input arguments
const std::string & (the parameter file name)
# return
-
```
> It reads in the parameter file name and initializes a ``Param`` object.
Copy and move assignment/semantics are disabled explicitly.

- **``Pipeline::assemble_grid``**
```
# input arguments
-
# return
-
```
> It is designed to initialize all grid instances with given parameter-set.


- **``Pipeline::assemble_breg``**
```
# input arguments
-
# return
-
```
> It is desinged to initialize the correct regular magnetic field model/type input.
It can also export the field if required.

- **``Pipeline::assemble_brnd``**
```
# input arguments
-
# return
-
```
> It is desinged to initialize the correct random magnetic field model/type input.
It can also export the field if required.

- **``Pipeline::assemble_tereg``**
```
# input arguments
-
# return
-
```
> It is desinged to initialize the correct regular thermal electron field model/type input.
It can also export the field if required.

- **``Pipeline::assemble_ternd``**
```
# input arguments
-
# return
-
```
> It is desinged to initialize the correct random thermal electron field model/type input.
It can also export the field if required.

- **``Pipeline::assemble_cre``**
```
# input arguments
-
# return
-
```
> It is desinged to initialize the correct cosmic-ray electron field model/type input.
It can also export the field if required.

- **``Pipeline::assemble_obs``**
```
# input arguments
-
# return
-
```
> It is desinged to initialize and calculate observables as required, with given fields.
It can also export the field if required.

## method list:

- **``main`` routine**

> The ``main`` routine is simple to read, so we directly attached the code below.
We emphasize that the ordering of calling ``Pipeline`` member functions are set on purpose.
The ``main`` function takes only one single argument, the input parameter-set file name.

```
int main(int /*argc*/, char **argv) {
#ifndef NTIMING
  auto tmr = std::make_unique<Timer>();
  tmr->start("main");
#endif
  const std::string input(argv[1]);
  auto run = std::make_unique<Pipeline>(input);
  run->assemble_grid();
  run->assemble_tereg();
  run->assemble_breg();
  run->assemble_ternd();
  run->assemble_brnd();
  run->assemble_cre();
  run->assemble_obs();
#ifndef NTIMING
  tmr->stop("main");
  tmr->print();
#endif
  return EXIT_SUCCESS;
}
```
