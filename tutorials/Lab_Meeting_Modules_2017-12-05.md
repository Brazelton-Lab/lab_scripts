
# Modules Tutorial

By: Alex Hyer (theonehyer@gmail.com)

This tutorial will cover how to use the `module` program on the Brazelton Lab Cluster. To follow along, you must login to said cluster using the crednetials provided to you by a Systems Administrator. Additionally, this tutorial will briefly cover BASH variables and environmental variables (EV) as understanding them is necessary for understanding `modules`.

## BASH Variables

BASH stands for Bourne Again SHell and is the most common terminal interface in \*nix systems. BASH is both a shell and a scripting language, meaning you can program in the tutorial. The following example makes three empty files using `touch`, uses a `for loop` to print files names starting with "test", and then deletes them with `rm`. *Note: a `for loop` is a programming construct. Thus, thus the following example uses code in the terminal.*


```bash
touch test1 test2 test3            # Create three empty files
for f in test*; do echo $f; done;  # Print file names, bad version of `ls`
rm test1 test2 test3               # Delete the files
```

    test1
    test2
    test3


Just like in a programming language, you can save and use variables in BASH. Variables can be assigned using the format `variableName=content`. Variables can be accessed using `$variableName` or `${variableName}`. The latter method is prefered as it allows you to insert variables into other strings/words. The following code block demonstrates these concepts: 


```bash
HELLO="WORLD"  # Define new variable HELLO
echo $HELLO    # Prints "WORLD", this is the most common example online
```

    WORLD



```bash
HELLO="WORLD"
echo ${HELLO}  # Prints "WORLD", this is the better habit
```

    WORLD



```bash
HELLO="WORLD"
echo $HELLOWORLD  # Prints nothing as there is no variable HELLOWORLD
```

    



```bash
HELLO="WORLD"
echo ${HELLO}WORLD  # Prints "WORLDWORLD"
```

    WORLDWORLD


Notice in the last two examples that `echo $HELLOWORLD` tries to print a variable called `HELLOWORLD` and `echo ${HELLO}WORLD` prints out the variable `HELLO` and adds "WORLD" to the end of it. This is why accessing variables via curly brackets is prefered.

BASH variables in the terminal can be used to save things like number of CPUs to use in programs, where to store log files, etc.

## Environmental Variables

Environmental variables are variables available to an entire operating system that affect how a system runs. For example, the `PATH` EV defines what folders contain executables (read programs) on the system and `HOME` points towards a users home directory.


```bash
echo $PATH
```

    /usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin



```bash
echo $HOME
```

    /home/ahyer


The `printenv` command will print all current environmental variables. There are quite a few EVs so this tends to produce a lot of output.

## Modules

The program `module` allows users and administrators to define EVs in file that can be loaded and unloaded for specific programs, projects, etc. Since EVs affect how a system operates, `module` allows users to define how they want things to work in their session and is thus highly flexible. The following sections cover two ways we use `module` on the Brazelton Lab Cluster. Before discussing these sepcific case, I will cover a few basics of `module`. *Note: this is not a comprehensive tutorial.*

### Basics

Modules are files that can be utilized by `module`. "Loading" a module will change your EVs by adding, deleting, or modifying them. To load a module, use (note that you can tab autocomplete module names):


```bash
module load <module name>
```

And unload them with:


```bash
module unload <module name>
```

Modules can also conflict with each other and `module` will refuse to load conflicting modules. Make sure to use ` module switch` or unload the one module before loading another if two modules conflict. A few more commands with brief explainations follow.


```bash
module avail       # Prints all available modules and where they're located
module list        # List currently loaded modules
module unload all  # Remove all loaded modules
module switch module1 module2  # Unload module1 and load module2
```

### Programs

Ever since 2017-05-26, all newly installed or updated programs on the Brazelton Lab Software and versioned under a `module`-based system. Under this system, every version of a program we install is kept and made accessible. For example, as of this writing, the cluster has `anvio` versions 2.0.2, 2.2.2, 2.3.2, and 2.4.0. By default, the cluster always uses the latest version of a program. However, older versions can be loaded via:


```bash
module load programs/<program name>-<version number>
```

Type the following commands on the cluster to see the system in action:


```bash
module list               # You likely have nothing loaded
srun anvi-interactive -v  # This prints anvi'o's latest version
module load programs/anvio-2.3.2
module list               # Notice the loaded module
srun anvi-interactive -v  # Note the version is 2.3.2
module switch programs/anvio-2.3.2 programs/anvio-2.2.2
module list
srun anvi-interactive -v  # Note the version is 2.2.2
module unload programs/anvio-2.2.2
module list               # You should have nothing loaded
srun anvi-interactive -v  # Note the version is now the latest version
```

Note that we're using `srun` here, which means `anvi-interactive` is running on a node and not `winogradsky`, yet `module` still works. SLURM sends your current EVs to the node when running a job. To use this system with `sbatch`, just place the `module load` command immediately after the `sbatch` options in an `sbatch` script. For example:


```bash
#! /bin/sh
#SBATCH --output test.output

module load programs/anvio-2.3.2
anvi-interactive -v
```

## mkproject

The program `mkproject` is a small but flexible BASH program I wrote that mimics `mkdir` except it creates an entire project's structure and generates a module and basic documentation for said project. The rest of this tutorial will focus on `mkproject` and will be structured so the reader can follow along line-by-line.

### Sweet and Simple

In it's simplest form, `mkproject` just requires a project name: `mkproject <project name>`. This creates a folder for the project, copies all files from the `/etc/mkproject/skel`, creates a module using `/etc/mkproject/TEMPLATE.module`, and edits the `README.md` file copied from `/etc/mkproject/skel` to generate semi-custom documentation. We will explore all these functions using the following:


```bash
mkdir tutorial        # A temporary directory for this tutorial
cd tutorial
mkproject basic       # Make a project named "basic"
ls                    # Note the firectory called "basic"
ls basic              # We can now see the structure of the project
ls basic/results      # Results has some additional directories
less basic/README.md  # This markdown file details the structure
module load basic     # Load the custom module for the project "basic"
                      # We'll print some EVs for "basic"
echo ${project}       # Points to project directory
echo ${data}          # Points to "data" folder in ${project}
echo ${results}       # Points to "results" folder in ${project}
printenv | grep -i -e '^project'  # These vars are in fact EVs
module list
ls ~/privatemodules   # Where modules are located, we'll open them later
module unload basic
```

There is an EV for each folder created automatically by `mkproject` as well as a `cpus` EV set to two. The project directory structure is derived from a pseudostandard used by many programmers when working on projects. While very few people reading this will be doing such work, the structure is well known and useful for consistency between projects. Since the project structure is detailed in `README.md`, I will not do so here.

If you don't want to strictly use this structure, you don't have to! In addition to simply editing the folders and module for your project, `mkproject` includes several builtin functions for tailoring things to your  style.

### .mkproject

If you have a directory called `${HOME}/.mkproject`, that directory is automatically used in place of `/etc/mkproject`. Note that `.mkproject` must meet the following requirements:

1. It must contain a file called `TEMPLATE.module`, this file can be empty.
2. It must contain a directory called `skel`, everything under `skel` will be copied to your project.
3. `skel` must have a file called `README.md`, this file can be empty.

This features allows you to overwrite the default project structure and module template in case you want your projects to always meet certain criteria. The following block will demonstrate how to make and use `.mkproject`.


```bash
mkproject -h                   # View the help for mkproject
demo="${HOME}/.mkproject"      # Hah, BASH variables! We come full circle.
mkdir -p ${demo}/skel          # -p makes parents along the way
touch ${demo}/TEMPLATE.module  # Required, can be empty
touch ${demo}/skel/README.md   # Required, can be empty
mkdir ${demo}/skel/test        # Example directory
mkproject config               # Make new project using the config we made
ls config
less config/README.md          # Note it is empty
module load config             # Note it breaks, let's find out why
echo ${project}  # Empty because we haven't properly configured a template
cat ~/privatemodules/config    # Note it is empty, hence we can't load it
```

The preceeding example is highly minimalistic, you'll want to make your own module template, README documentation, and add a more robust directory structure if you choose to use `.mkproject`. I recommend copying `TEMPALTE.module` and `README.md` to your `.mkproject` and editing from there. Both files have tutorials on hwo to use them in the file and already contain everything you'll need, you just have to edit them. Let's do that now:


```bash
less /etc/mkprokect/skell/README.md
less /etc/mkproject/TEMPLATE.module
less ~/privatemodules/basic          # Let's look at one we've made
```

You may need to play around with your these files to get the hang of them, but it'll be worth your while.

### Using the Default

Even if you make `.mkproject` to create a custom project structure, you may still want to use the defaults stucture from `/etc/mkproject` on occasion. To do so, simply use the flag `-d` when running `mkproject`. Let's try it:


```bash
mkproject -d def  # New project using original template, default is reserved
ls def            # See how directory structure is the same as "basic"
less def/README.md
module load def
echo ${project}
module unload def
```

### Multiple Custom Project Structures

While `.mkproject` overrides the default project structure to suite your tastes, you may want multiple different structures for different types of projects, e.g. 16S or metagenomics. The `-c` options lets you specify any directory you want to be the skeleton for your project, however, it must follow the same rules as `.mkproject`. The following examples illustrates how to make custom configurations outside of `.mkproject`:


```bash
newconfig="${HOME}/.test"          # This will be the custom config
mkdir -p ${newconfig}/skel
touch ${newconfig}/TEMPLATE.module
touch ${newconfig}/skel/README.md
mkdir ${newconfig}/skel/test2      # Example directory, different from above
mkproject -c ~/.test custom        # New project using custom config
ls custom
```

And that's `mkproject`!

## Final Notes

1. Modules can "stack" as long as they don't conflict. You could have a module for your project and additional modules for different groups of sample or analyses.
2. `module` is very flexible, be creative!
3. You can save program versions per project by using the lines of code from the programs modules. Simply copy the lines at the bottom of the file into your project's module.
4. I may add a "freeze environment" option to `mkproject` to automate the preceeding point. No promises.
5. If you make a project with the same name as a previous project, the newer project will not get a module.

## Cleanup

Let's delete all the files and folders we made in this tutorial.


```bash
cd ../
rm -rf tutorial .mkproject .test
rm ~/privatemodules/basic ~/privatemodules/def ~/privatemodules/config
```
