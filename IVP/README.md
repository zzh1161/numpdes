# Project 3 文件与程序说明

## 文件结构说明
+ `bin` 目录下是测试和绘图生成的可执行文件
+ `include` 目录下是所有的头文件
  + `Info.hpp` 里声明了精度要求的宏定义
  + `Vec.h` 定义了向量，所有的计算都基于此
  + `IncAllMethod.h` 包含了所有的方法的头文件，便于测试程序引用
  + `GenericFactory.hpp` 里定义了泛型抽象单体工厂基类
  + `TimeIntegrator.hpp` 是所有方法的抽象基类
  + `json.hpp` 是一个header-only的json解释包，地址为https://github.com/nlohmann/json
  + 其余头文件均为相应方法类的声明和实现
+ `output` 目录下是程序的输出
  + `test*.txt` 这些文件是相应的测试文件输出的误差、收敛率、运行时间等
  + `plot.txt` 是用于matlab绘图的格式化输出
+ `report` 目录下是报告文档
+ `src` 目录下是使用对象工厂模式调用作业要求中的三种方法绘图，并分析运行时间的程序
  + `matplot` 目录下是绘图的matlab程序，绘制的图像在`report/image` 下
  + `inputfile` 下是json格式的输入文件
+ `test` 目录下是所有的测试程序，分别测试了作业要求中的8种方法
  + `inputfile` 下是json格式的输入文件
  + `catch.hpp` 是一个header-only的C++单元测试框架，地址为https://github.com/catchorg/Catch2

## 编译与运行

代码作者本地具有以下环境和软件：

+ Ubuntu 20.04
+ g++(GCC) 11.2.0 (C++20)
+ GNU Make 4.1
+ TeX Live 2021

编译：

+ 输入`make` 编译`test`和`src`目录下所有的cpp文件
+ 输入`make story` 编译报告文档

运行：

+ 输入`make run` 将执行结果输出到指定文件中
+ 输入`make clean` 清除所有可执行文件和输出的txt文件

## 程序报告

报告见[这里](./report/report.pdf)