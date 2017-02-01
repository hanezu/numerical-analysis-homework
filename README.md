# 2016 Numerical Analysis Homework

Created 004 for problem 1 and 005 for problem 2.
003 was also created for a simpler version of problem 1, but it was not completed.

The notes and report questions are also included under the root directory.

# To run a particular program 

Execute run.py under python3 to run a particular program by calling its folder name.

ALL python scripts are written in python3, and may not be able to work in python2.

## options:

  -m: make the c file.


  -f: take the next argument as the folder of code to run.

      enter the full name, or the first 3 digits would be enough.

      Therefore the folder name (and the c file name) must obey the naming rule.

  -p: plot the result picture when the program terminate.  


# example of using run.py

python3 run.py -m -f 004 -p

python3 run.py -m -f 005

# To change the sweeping parameters

To change the parameter of the number of splitting in space (h) or in time (tau), change the setting in the c file. These setting are mainly included in the function "main".

It is also possible to run the program with default setting of parameter zone and splitting.

Below is the original README file.

# 京都大学 2016 年度後期「数値解析」 サンプルプログラム集
- 担当教員：Karel Svadlenka
- プログラム作成者：宇田 智紀

このリポジトリでは，数値解析の授業で紹介された偏微分方程式の差分解法のサンプルプログラムを公開しています．

すべてのプログラムは GNU Make および GCC にてコンパイルできます．また，可視化は gnuplot を使うのが便利です．

```
$ make
$ ./program.out > output.dat
$ gnuplot -e "splot 'output.dat'"
```

パブリックドメインですので，自由に改変・再配布等して試してもらってかまいません．ライセンスの詳細はLICENSE ファイルを参照してください．

