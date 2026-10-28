# 从公式推导到同精度 C++ 代码

完整的 JSON 输入、源代码 diff 输出和英文复现说明见 [README.en.md](README.en.md)。

从仓库根目录一条命令生成可应用的补丁：

```console
python tools/ofp_codegen/repair.py --spec tools/ofp_codegen/examples/cos_x2.input.json --project-root . --output runs/cos_x2-patch
```

`repair.py` 会读取输入指定的源函数，输出 `patch.diff`、预览代码、生成的头文件、
系数记录和 `manifest.json`，不会直接修改原函数。输入中的展开方案和代码定位由使用者提供。
`examples/project_cos_x2.input.json` 可以针对项目的 `inputs/targets.cpp` 生成补丁。
这是新增的自动化工具，不表示论文已有的 57 个补丁原先均由该程序自动产生。

这个目录可以整体复制到另一台电脑独立运行。核心生成器 `generate.py` 不读取 OFP-Repair 的
`inputs/targets.cpp`、`inputs/target.h`、已有 patch 或 57 例的审计结果。
它接收一个数学表达式和明确的展开设置，自动推导系数，生成同类型 C++ 函数，
同时保存推导记录。目前它位于 OFP-Repair 仓库内，**不是独立的 Git 仓库**。

## 输入与输出

必须提供：数学表达式、展开点、最高次数、适用区间、函数名及输出路径。
变量名默认 `x`，浮点类型默认 `double`，也支持 `float`。

自动执行的步骤：

1. 用受限数学语法读取表达式，数值字面量按精确有理数解释。
2. 把展开点舍入到目标浮点类型，并在该浮点数的**精确数学值**处展开。
3. 在截断前检查局部亚纯性，再计算左右展开并提取多项式系数。
4. 把系数舍入到目标浮点类型，生成 Horner 形式的 C++ 代码。
5. 生成区间判断：区间内调用多项式，区间外调用用户传入的原函数。
6. 保存 `.hpp` 和同名 `.json`，记录符号系数、浮点常量、位模式、版本和 SHA-256。

输出 C++ 只使用所选的 `float` 或 `double` 运算，不使用运行时任意精度库。
离线符号推导和常量计算可以使用高精度；这不等于提升生成函数的浮点类型。
提高精度来自更稳定的表达式，能否达到目标误差仍须验证。

## 快速运行

需要 Python 3.10 或更新版本。进入此目录后安装依赖：

```console
python -m pip install -r requirements.txt
```

生成 `(1-cos(x))/x**2` 在零点附近的 10 次多项式：

```console
python generate.py --expression "(1-cos(x))/x**2" --center 0 --degree 10 --interval -0.01 0.01 --dtype double --name cos_x2 --output examples/cos_x2.hpp
```

符号结果为：

```text
1/2 - u**2/24 + u**4/720 - u**6/40320 + u**8/3628800 - u**10/479001600
u = x
```

`examples/demo.cpp` 展示如何把生成代码接到原函数上：

```cpp
#include <cmath>
#include "cos_x2.hpp"

double original(double x) {
    return (1.0 - std::cos(x)) / (x * x);
}

double repaired(double x) {
    return ofp_codegen::cos_x2_repaired(x, &original);
}
```

编译示例（Clang/GCC，C++17）：

```console
c++ -std=c++17 -O2 -fno-fast-math -ffp-contract=off examples/demo.cpp -o demo
```

运行生成器不需要编译器；编译和数值验证才需要编译器。Docker 不是程序依赖。
生成 `float` 版本只需设置 `--dtype float`，其常量带 `f` 后缀，接口也使用 `float`。

## 怎样复查

```console
python validate.py --cxx clang++ --output validation-output
```

`--cxx` 可以是 Clang/GCC 编译器路径；验证脚本的命令行参数不兼容 MSVC `cl`。
输出目录必须不存在或为空，防止覆盖此前的验证证据。

验证脚本实际调用生成器、编译并运行生成的 C++，用 mpmath 的 100 位参考值核对
**实际 C++ 输入和输出**，而不是只在 Python 中代入同一多项式。

- 三组表达式：`(1-cos(x))/x**2`、`(x-1)/(x**2-1)`、在非零中心展开的 `sin(x)`。
- 每组分别生成 `double` 和 `float`，共六个实现；共 54 个区间内数值检查、30 个回退检查。
- 覆盖中心、端点、相邻浮点数、NaN 和无穷；静态检查返回类型，严格编译检查隐式精度提升。
- 拒绝极点、分支点、不可解析表达式及非法语法，包含低阶截断可能掩盖的分支项。
- 检查大小写不同的函数名可以在同一 C++ 文件中同时包含。
- 检查原始消减表达式与生成实现的实际误差改善。

每个步骤保留 `command.json`、`result.json`、`stdout.log` 和 `stderr.log`；
总结果在 `run-report.json`，包含逐点参考值、绝对误差、相对误差和检查状态。
生成的头文件、符号推导 JSON、C++ 测试程序和可执行文件也保留在输出目录。

## 适用范围

这是一个**给定展开方案后的自动推导与代码生成器**。以下设置仍由调用者提供：
数学表达式、展开变量、中心、阶数及适用区间。程序不自动决定它们，也不自动识别任意
C++ 函数的控制流、库调用语义或失稳位置。目前没有自动批量导出原项目 57 例完整 patch 的入口。

语法支持 `+ - * / **`、常量 `pi` 和 `E`，以及单参数函数
`sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh exp log sqrt erf erfc Abs`。
语法识别不代表每个函数在每个点都能生成：SymPy 必须先明确判定局部亚纯性，
不确定就拒绝。当前 `erf`、`erfc` 等解析函数也可能因判定能力不足被保守拒绝。
负幂、分数幂或左右不一致的截断结果会被拒绝。当前最高次数为 64。

需要倒数换元、对数/平方根外包、特殊函数映射、尾部渐近展开或 Newton 迭代的案例，
还需要增加对应的策略模块。已有 57 例审计中的这些人工设置不会被本工具自动推断出来。

## 数值语义与限制

- 区间端点向内舍入，比较是闭区间比较。NaN、正负无穷及区间外输入调用原函数；
  调用者应传入非空、同输入输出类型的原函数指针。
- 多项式在可去奇点处返回连续延拓值，例如本例零点返回 `0.5`。
  不保证保留原函数的 errno、浮点异常或有符号零行为。
- 本工具依赖 SymPy 简化后的数学表达式，不证明原始表达式定义域完全保持，
  也不证明所给整个区间都解析或达到误差要求。区间不能仅凭代码生成成功就直接采用。
- 有理系数用整数算法做最近偶数舍入；其他常量在 120 和 180 位计算后须得到相同
  浮点数，否则拒绝。这是交叉精度检查，不是超越常量正确舍入的形式证明。
- 系数下溢会记录在 JSON 中；生成成功不保证多项式足够准确。
- `.json` 中的 `generated_unvalidated_interval` 表示生成成功、区间尚未获精度认证。
  有限采样和编译检查也不构成整个区间的误差上界证明。
- 验证使用 IEEE binary32/binary64 环境，禁用 fast-math 和 FMA 收缩。
  改变编译器浮点优化或执行环境后应重新验证。

## 文件

| 文件 | 用途 |
| --- | --- |
| `generate.py` | 数学输入 → 符号展开 → `.hpp` 与推导 `.json` |
| `repair.py` | JSON 规范 → 系数推导、源代码 diff、预览和 manifest |
| `pipeline_checks.py` | LF/CRLF diff 应用、编译和误差改善的端到端验证 |
| `validate.py` | 实际 C++ 编译、数值与回退验证，并保留全过程 |
| `requirements.txt` | 固定版本的 Python 依赖 |
| `examples/cos_x2.hpp`、`.json` | double 示例与符号推导记录 |
| `examples/sin_float.hpp`、`.json` | float 非零中心展开示例 |
| `examples/demo.cpp` | 用原函数回调接入生成代码的示例 |
