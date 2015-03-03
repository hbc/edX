### Introduction to Unix [SCREENCAST]

####Introduction and the terminal

A short slide deck as part of the screencast to introduce UNIX/Linux and command line concepts

####Finding your way around the file system
>Note – All UNIX commands are case sensitive

<pre>
$ ls      # list all the files and subdirectories in the current directory

$ ls unix_exercise/       # list all the files and subdirectories in the directory unix_exercise

$ ls unix_exercise/readme.txt       # list the file readme.txt within unix_exercise
</pre>

The colors in the listing of the items in the current directory indicate what type of file/folder it is.
<pre>
$ ls –F unix_exercise/         # “-F” is an argument for the command ls, we will talk more about arguments momentarily
</pre>

In the listing -

* All entries with a `/` at the end <-> directories
* All entries with a `*` at the end <-> executables
* All entries with a `@` at the end <-> linked directory or file (“symbolic link” or a “sym link”). This is a directory/file linked to your home page, but is actually a branch elsewhere on the directory tree. This is useful for easy access.
* The rest of the entries are files.

Lets make a new file
<pre>
$ touch testfile          # touch lets you create a new, empty file

$ ls
</pre>

You can see that you have created a simple file, now lets remove or delete that file we just created

<pre>
$ rm testfile        # rm = remove file. 
</pre>
>
Note – File naming convention in UNIX has certain features:
>
>* Use only letters (upper- and lower-case), numbers from 0 to 9, a dot (.), underscore (_), hyphen (-). Do not use “space” in file names.
>
>* Avoid other characters, as they may have special meaning to either Linux, or to the application you are trying to run.
>
>* Extensions are commonly used to denote the type of file, but are not necessary. On the command line you always explicitly specify a program that is supposed to open or do something specific to the file.
>
>* The dot (.) does not have any special meaning in Linux file names. (Remember that the dot does have a special meaning in the context of the command line though.)


What is your location right now in the context of the directory structure? You are in your home directory… But, “where” in the tree of the directory structure is your home directory?

<pre>
$ pwd          # pwd = print working directory
</pre>
Lets change our working directory to unix_exercise
<pre>
$ cd unix_exercise        # cd = change directory
</pre>

A *Relative* path, is a path to the location of your file of interest relative to your working directory (or location) (e.g. `../file.txt`), whereas the *Full* path is the location of the file starting with the root directory (e.g. `/home/vagrant/text_files/file.txt`). 

In the `cd` command above you used a relative path.
<pre>
$ cd      # no matter where you are, if you say just “cd”, the OS returns you back to your home directory

$ cd /home/vagrant/unix_exercise

$ pwd

$ cd

$ cd unix_exercise

$ pwd
</pre>

#### Manipulating files and directories
Making new directories is very simple, lets make a new one called new_dir.
<pre>
$ mkdir new_dir      # mkdir = make new directory
</pre>
How is the following command different from the previous one?
<pre>
$ mkdir new dir      # two more directories new and dir will be created because of naming conventions
</pre>

<pre>
$ rm new             # “rm: cannot remove `new`: Is a directory”
</pre>

We need an argument to go along with rm to enable it to delete directories. Man pages (manual for each command) are very helpful in figuring out the arguments you can use with specific commands. (Other than man pages, the internet is a good resource for more information about commands, but be discerning.) Arguments help you give the command special instructions and do more specific tasks.

<pre>
$ man rm             # man = manual for a specific UNIX command. typing the letter q (lower case) will get you back to the command line

$ rm –r new

$ man ls
</pre>

Let's backup the unix_exercise directory; we can copy the whole directory to a new directory:
<pre> 
$ cp –r unix_exercise unix_exercise_backup         # cp = copy file or directory (-r). The first directory or file name you state after the command is what is being copied; the second file name is what you are copying to. When you use copy, if the a directory/file, that you state second, doesn’t already exist it is created
</pre>

<pre>
$ cd unix_exercise/
</pre>

Create a new file using touch, move it to the home directory and then rename it:
<pre>
$ touch new_file.txt

$ mv new_file.txt /home/vagrant/     # moved.

$ cd /home/vagrant/

$ mv new_file.txt home_new_file.txt      # renamed! mv can move and rename files.
</pre>

>
Note – As we start learning more about manipulating files and directories, one important thing to keep in mind is that unlike Windows and Mac OS, this OS will not check with you before replacing a file. E.g. if you already have a file named foo.txt, and you give the command `cp boo.txt foo.txt`, all your original information in foo.txt will be lost.

#### Examining file contents

So far we have learned to move files around, and do basic file and directory manipulations. Next we’ll learn about how to look at the content of a file.

The commands, `cat`, `head` and `tail` will print the contents of the file onto the monitor. `cat` will print ALL the contents of a file onto your terminal window, so be aware of this for huge files. 

`head` and `tail` will show only the number of lines you want to see (default is 10 lines):
<pre>
$ cat readme.txt          # cat = catenate, prints the whole file

$ cd sequence/

$ head chr4.fa       # prints on the screen the first 10 lines of the file

$ tail chr4.fa       # prints on the screen the last 10 lines of the file
</pre>

The commands `less` and `more` allow you to quickly take a look inside a file. With `less` you can use the arrow keys to go up and down, pressing the “q” key will get you back to the command prompt. This is similar to what you encountered with the `man` command. `more` is not as good for large files, since it loads the whole file into memory, and also it doesn’t let you go backwards within the file. So we will stick to using `less`.
<pre>
$ less chr4.fa
</pre>

#### Output redirection (>, >> and |)

What if we wanted to collect the top and bottom 50 genes from genelist1.txt in the genelists directory, and make a new file with the 100 genes. 
<pre>
$ cd

$ cd unix_exercise/genelists/

$ head –n 50 genelist1.txt > genelist_test1.txt          # “>” redirects output to specified file, but if the file already existed it overwrites the contents!

$ tail –n 50 genelist1.txt > genelist_test2.txt

$ cat genelist1_test1.txt genelist_test2.txt > genelist1_test_combined.txt           # the “cat” command will print the contents of both files in the order you list them, and in this case you have redirected the merged output from the 2 files to a new file instead of the terminal
</pre>

We used 3 steps above to get the combined file, instead we could have done it in 2 steps (see below).
<pre>
$ head –n 50 genelist1.txt > genelist1_test_combined.txt           # this will overwrite the file

$ tail –n 50 genelist1.txt >> genelist1_test_combined.txt          # “>>” redirects the output to specified file, however it appends the new content to the end of the file and does not overwrite it
</pre>

Pipes or `|` is a very handy UNIX tool to string together several commands into 1 command line. Basically, it takes the output of one command and “pipes” it into the next command as input.

What if we also needed to make sure that the new document had the data sorted alphabetically?
<pre>
$ cat genelist1_test1.txt genelist1_test2.txt | sort > genelist1_test_combined_sorted.txt       #sort = sorts data as you specify in the arguments, default is alphanumeric in ascending order

$ head genelist1_test_combined*			# the asterisk "*" is a wildcard and can be used in place of 1 or multiple characters
</pre>

#### Permissions
UNIX is a multiuser system, and to maintain privacy and security, most users can only access a small subset of all the files.

* You are the owner of every file and directory that is under your home directory.
* The system administrator (sys admin) or other users can determine what else you have access to.
* Each file and directory has associated “permissions” for different types of access; reading, writing and executing (scripts or programs).
* You are allowed to change the permissions of any file or directory you “own” and in some cases a file or a directory that you have access to as part of a “group” (co-ownership). 

<pre>
$ ls -l /home/vagrant/unix_exercise/
</pre>

<figure>
<a href="/Users/rkhetani/Git_Repos/edX/images/ScreenShot.png">
<img src="/Users/rkhetani/Git_Repos/edX/images/ScreenShot.png">
</a>
</figure>

`d`: directory (or `-` if file); 

`r`: read permission; 

`w`: write permission; 

`x`: execute permission (or permission to `cd` if it is a directory); 

`-`: no permission.       

`drwxr-xr--` can be divided into `d` `rwx` `r-x` `r--`, and it means the following:

* owner (u) has `rwx` read, write and execute permissions for the directory

* group (g) has `r-x` only read and execute permissions for the directory

* others (o) has `r--` only read permission for the directory

How do you set or change these permissions? 
<pre>
$ cd ../

$ chmod -R o-rwx sequence/     # others have no read, write or execute permissions for any this directory or any file within

$ ls -lh 

$ chmod u+rwx hello_world.sh

$ ls -lh

$ chmod -R 764 sequence/       # same as “chmod –R u+rwx,g+rw,o+r”

$ ls -lh
</pre>
A “sticky bit” is applied to shared directories to protect files such that only the owner has the ability to change permissions. 

`chown` and `chgrp` are commands that let you change owner and groups respectively, but you need to start out with correct permissions to be able to execute these on a file/directory.
