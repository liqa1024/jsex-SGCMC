# Requirement

- [jse](https://github.com/liqa1024/jse/) (latest)

- C & C++ Compiler

  For Windows, [MSVC](https://visualstudio.microsoft.com/vs/features/cplusplus/) is recommended

  For Linux, [GCC](https://gcc.gnu.org/) is recommended

- MPI Development Environment

  For Windows, [Microsoft MPI](https://www.microsoft.com/download/details.aspx?id=105289) is recommended
  (both `msmpisdk.msi` and `msmpisetup.exe` are required)

  For Linux like Ubuntu, you can use `sudo apt install libopenmpi-dev`


# Usage

1. Ensure that the environment requirements mentioned above are met and that `jse` has finished compiling the relevant JNI libraries:

   ```shell
   jse --jnibuild
   ```

2. Ensure that the `lmp` folder is in the `jse` Groovy classpath. You can choose one of the following methods:

    * `jse` adds the current working directory to the classpath, so you can simply ensure that the `lmp` folder exists in the current running directory.

    * Place the `lmp` folder inside the `lib/groovy` directory of `jse`. `jse` automatically adds this directory to the classpath:

      ```text
      ├─lib
      │ ├─groovy
      │ │ └─lib
      │ │   ├─FixSgcmcNNAP.groovy
      │ │   ├─FixFlipNNAP.groovy
      │ │   └─FixFlipNNAPScaled.groovy
      │ └─jse-all.jar
      └─jse
      ```

    * Add the directory containing `lmp` to the `JSE_GROOVY_EXLIB_DIRS` environment variable:

      ```shell
      export JSE_GROOVY_EXLIB_DIRS="path/to/your/package/dir:$JSE_GROOVY_EXLIB_DIRS"
      ```

3. Run LAMMPS via `jse -lmp`. For example:

   ```shell
   mpiexec -np 8 jse -lmp -in flipFe.lmpin
   ```

---

As a many-body interaction, NNAP requires double the cutoff radius to strictly calculate the energy difference caused by species changes. Therefore, you need to set a larger communications cutoff radius using a command similar to:

```text
comm_modify cutoff 14.0
```

($r_{c} = 6.0,~~6.0 \times 2.0 + 2.0 = 14.0$). A simpler method is to configure the fix and run it serially; the program will output a suggested setting prompt. For specific settings, please refer to the example file `flipFe.lmpin` and the comments in the source code.


# Citation

_TODO_


# License

Scripts on this repository are licensed under the **GNU GPL v3**.
See [LICENSE](LICENSE) for details.

