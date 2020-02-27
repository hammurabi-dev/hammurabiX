## hammurabi X Timer class

The ``Timer`` class is designed for recording executing time of an arbitrary part/fraction of the simulation routine.

- [header file](https://bitbucket.org/hammurabicode/hamx/src/master/include/timer.h)

## function list:

- **``Timer::start``**
```
# input arguments
std::string & (key word for timing)
# return
-
```
> It records the key word and start timing, until the ``stop`` function is called.

- **``Timer::stop``**
```
# input arguments
std::string & (key word for timing)
# return
-
```
> It stops timing of the corresponding key word

- **``Timer::print``**
```
# input argumnet
std::string & (key word for timing)
# return
-
```
> It prints the timing result of given key word, if the input argument is left empty, it by default prints all recorded key word timing results.

## method list:

- **timing with key word**

> In executing the routine, we may want to record the time for running a specific part of the whole pipeline.
We design the ``Timer`` class which in principle hosts a C++ standard library map structure, where the key word acts as a map entry, and the timing result is registered as content under each key word.
In this way, we can mix several timing actions in a single routine execution and read out the results conveniently.
