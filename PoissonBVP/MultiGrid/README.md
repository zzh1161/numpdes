# 文件与程序说明

本次程序完成了除附加项(c)以外的所有任务.

## 文件结构说明

+ `bin` 目录下是生成的可执行文件
+ `include` 目录下是程序所需头文件
  + `BeforeAll.hpp` 是前置信息声明
  + `BoundaryCondi.hpp` 是边值条件类声明和定义
  + `Grid.hpp` 声明了Grid类，通过将其特化实现不同维度的多重网格方法
  + `GridDimOne.hpp` 声明并定义了一维网格，并实现了相应的多重网格方法
  + `GridDimTwo.hpp` 声明并定义了二维网格，并实现了相应的多重网格方法
  + `json.hpp` 是一个header-only的json解释包，地址为https://github.com/nlohmann/json
+ `output` 目录下是程序的输出
  + `matplot` 目录下是绘图所用的matlab程序
  + `figure` 是绘制的图片等
+ `report` 目录下是本次作业的报告
  + `reportAB.pdf` 是本次作业(A)(B)部分的报告
  + `theory2D.pdf` 是附加项(b)中的二维多重网格的理论分析部分
+ `src` 目录下是主程序的cpp文件
  + `main1D.cpp` 是(A)部分的程序
  + `main2D.cpp` 是(B)部分的程序
  + `display.cpp` 是格式化输出以便绘图的程序
+ `test` 目录下是测试文件
  + `catch.hpp` 是一个header-only的C++单元测试框架，地址为https://github.com/catchorg/Catch2

## 编译与运行

代码作者的本地具有以下环境和软件：

+ Ubuntu 18.04
+ g++(GCC) 11.2.0 (C++20)
+ Lapack 3.10.0
+ GNU Make 4.1

编译：

+ 输入`make` 生成主程序部分的可执行文件
+ 输入`make test` 生成测试
+ 输入`make clean` 清除生成的可执行文件与输出的文本文件

运行：

+ 输入`make run` 将可执行程序的输出定向到指定文件中

## 设计说明

### 类设计

+ 首先设计了Grid类，通过类特化实现不同维度的Grid类及相应的多重网格方法

  ```cpp
  template<int dim>
  class Grid{};
  ```

+ 一维的Grid类里有成员变量一维数组points\_用于表示选取的一维网格上的点，一维数组RightHandSide\_用于表示求解的线性方程组的右端项.

+ 二维Grid类也具有同一维相似的成员变量，只不过替换为了二维数组.

### 实现思路

+ 首先，Dirichlet、Neumann、Mixed三种边值条件其实可归结为一种，即Mixed. 通过获取第三边值条件中的系数a、b，

  $$
  a\frac{\partial U}{\partial \vec{n}}+bU = \sigma
  $$
  
  利用边界点附近的三个点推导二阶准确表达式，

  $$
  \left(\frac{3a}{2h}+b \right)U_0 - \frac{2a}{h}U_1 + \frac{a}{2h}U_2 = \sigma
  $$
  

+ 函数`V_cycle()` 和`FMG()` 设定为类的私有成员函数，外部调用通过函数`Poisson_BVP_Multigrid()` 预留接口，选择方法实现

+ 函数`analytic_error()` 表示解与解析解的误差的无穷范数大小

  函数`discrete_error()` 通过直接解法求得离散后的准确解，并计算解与离散准确解的误差的无穷范数

## 程序报告

(A)(B)部分的分析报告请见[./report/reportAB.pdf](./report/reportAB.pdf)

附加项(b)部分的理论分析请见[./report/theory2D.pdf](./report/theory2D.pdf)

