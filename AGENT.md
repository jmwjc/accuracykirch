# AGENT.md — 本仓库的撰写格式约定

本文件约束 `manuscript/` 下所有 `.tex` 的写法，修订与新增内容一律遵守。

## 硬性规则

- 正文**一句话一行**：每句话单独占一行，句末换行。
- 不得在一句话中间硬换行，也不得把两句话挤在同一行。
- 公式环境用缩进表达层次：`\begin{equation}` 顶格，环境内每一行缩进一个 tab。
- 嵌套环境（如 `cases`、`aligned`）各多缩进一个 tab，`\end{...}` 与其 `\begin{...}` 对齐。
- 公式过宽时拆成多行（`aligned` / `split` / `gathered`），不要靠缩小字号或挤压间距硬塞。
- 括号默认用圆括号 `\left( ... \right)`，按需嵌套；方括号 `[]` **仅**用于表示跳跃或角点两侧之差，不作普通分组。
- 跳跃量用 `\left[ ... \right]` 表示，如 $[w_{,ij}n_is_j]$；不用上标 `^{+}` / `^{-}` 写差值。
- 不用 `\nabla`：微分一律用逗号下标表示，如 $\nabla^4w\to w_{,iijj}$、$(\nabla^2w)_{,i}\to w_{,jji}$。
- 公式里相邻变量之间**加空格**，如 `s_i s_j`、`n_i n_j`、`w_{,ij} n_i n_j`、`\Psi_I(x) d_I`。
- **行间公式末尾不加标点**：句号、逗号一律不写（`\end{equation}` 前的最后一行不以 `.` 或 `,` 结尾）。
- 空格直接写在数学模式里，不要把变量间距靠 `\,` 或 `\;` 硬拼。
- 加粗一律用 `\boldsymbol{...}`，不用 `\bm{}`、`\mathbf{}`。

## 示例

正确：

```latex
Continuity follows directly from the Cauchy--Schwarz inequality.
Coercivity follows from the positive definiteness of the elasticity tensor together with the Poincar\'e-type inequality.
\begin{equation}\label{continuity}
	\left| a(\delta w, w) \right| \le \gamma \Vert \delta w \Vert_{H_2} \Vert w \Vert_{H_2}, \quad \forall \delta w, w \in H_2
\end{equation}
```

错误：

```latex
Continuity follows from the Cauchy--Schwarz inequality. Coercivity follows from the positive definiteness of the elasticity tensor.
\begin{equation}\label{continuity}
\left| a(\delta w, w) \right| \le \gamma \Vert \delta w \Vert_{H_2} \Vert w \Vert_{H_2}, \quad \forall \delta w, w \in H_2
\end{equation}
```

## 作者写作风格

下列风格取自作者亲自改定的稿件，**优先于任何通用惯例**。

- 定义新符号时给**独立编号公式**，不合并进 `gathered`、不写在正文行内。
- 引定义的固定句式：`where $X$ denotes ..., and it can be expressed by A, B and C as:`，紧跟公式，公式后接 `and $Y$ ...`。
- 多个小公式宁拆成多条 `equation`，不堆进一个 `aligned`。
- 性质与条件用 `itemize` + `\item\em{...}` 列出，每个 `\item` 后接编号公式，列表后可补一句收束。
- `H_2` 直接按 Sobolev 空间 $H^{2}(\Omega)$ 使用，不解释其是否含本质边界条件；本质边界条件后续由弱形式施加（避免刚度矩阵奇异），此处不作强调。
- 段落开头惯用 `where ...`、`It is noted that ...`、`As a result, ...`、`In this context, ...`、`However, ...`。
- 角点集写 `c_P`；角点在文字中称 `corners of boundary curve $\Gamma^{c}$`。
- 带上划线的符号表示给定值或材料常数（`\bar w`、`\bar V`、`\bar D`），带括号的函数名表示由 $w$ 导出的量（`\theta_n(w)`、`V(w)`、`M_n(w)`、`P(w)`）。
- **不改写作者已有句式和用词**，只做点对点替换；发现疑似笔误先提出，不擅自纠正（例如目前保留 `Possion` 写法）。
- **不得删除作者新增的公式、段落或文件内容**：只做记号统一、错字修正与必要的上下文补写；确需改动时先提出，由作者决定。

## 编译（唯一事实来源：Helix 的 texlab 配置）

编译方式由 `~/.config/helix/languages.toml` 中 `[language-server.texlab.config.texlab.build]` 定义，等价命令为：

```bash
cd manuscript && latexmk -xelatex -auxdir=.aux -interaction=nonstopmode -file-line-error -synctex=1 main.tex
```

- 引擎固定 XeLaTeX，参数固定为上面这一组；不得再加入 `-pdf`、`-draftmode`、`-outdir` 等开关（`-pdf` 与 `-xelatex` 冲突、`-outdir` 会让 texlab 找不到文件）。阅读器相关的开关（sioyek 的 `--nofocus`、galleypdf 的 `-g`）**只能**出现在 `forwardSearch.args` 里：写进 `build.args` 会被当成 latexmk 参数，报 `Bad options specified`，整次编译直接失败。
- 前向搜索固定用 galleypdf：`executable = "galleypdf"`，`args = ["forward", "-g", "-l", "%l", "-s", "%f", "%p"]`，其中 `-g` 表示不把 Galley 提到前台。**不要加 `%c`**：texlab 对 `%c` 的替换会把其后的参数粘成一个（实测收到 `20--ccc--s//…`），列号也非必需。占位符只用 `%f`、`%l`、`%p`。
- 中间文件固定落 `manuscript/.aux/`（`auxDirectory = ".aux"`）；`main.pdf` 与 `main.synctex.gz` 固定落 `manuscript/`（`pdfDirectory = "."`）。这两处位置一旦改动，正向/反向搜索与 texlab 诊断即失效。
- 日常校查用同一套 XeLaTeX 参数、加 `-no-pdf` 跑两遍，**不碰 `main.pdf`**：

```bash
cd manuscript && xelatex -no-pdf -output-directory=.aux -interaction=nonstopmode -file-line-error -synctex=1 main.tex && xelatex -no-pdf -output-directory=.aux -interaction=nonstopmode -file-line-error -synctex=1 main.tex
```

- 注意：`-auxdir=.aux` 是 **latexmk** 的选项，`xelatex` 本身不认（会报 `unrecognized option`）；直调引擎时必须用 `-output-directory=.aux`。

- 需要出 PDF（Helix 保存时 `onSave = true` 会自动触发，或人工核对版式）才用上面的 `latexmk` 命令。
- `main.pdf` 已纳入版本库：Helix 的 onSave 重编会把它改脏，提交前执行 `git restore manuscript/main.pdf`，除非本次确实要更新 PDF。
- 校查判据只读 `manuscript/.aux/main.log`：`grep -n "undefined\|Undefined" .aux/main.log` 必须无输出，`grep -c Overfull .aux/main.log` 不得比上次多。
- 本节与 `languages.toml` 必须逐字对应：改动任一侧（引擎、参数、`auxDirectory`、`pdfDirectory`）都要同步另一侧。
- 报错定位以上面的 `-file-line-error` 输出为准（`文件:行: 错误`），与 Helix 里 texlab 报诊断的行号一致。

## 其他约定

- 一节一个文件，置于 `manuscript/sections/`，由 `manuscript/main.tex` 用 `\input` 引入。
- label 命名用前缀区分：`sec:`、`fig:`、`tab:`；公式 label 用**语义名**（如 `weakform`、`nitsche_coercivity`、`defect_bound`），不用 `eq:数字` 这类编号。
- 交叉引用统一用 `\eqref{}` 引用公式、`\ref{}` 引用章节与图表。
- 参考文献为手写 `thebibliography`（`manuscript/sections/references.tex`），新增文献须手动补 `\bibitem`；`references.bib` 仅作资料留存，不参与编译。
- 公式改动后会整体重编号，收尾时须回查全篇 `\eqref`。
- 编译一律按「编译（唯一事实来源：Helix 的 texlab 配置）」一节执行，不再使用 `pdflatex`。
- 每次改完确认无 `Undefined` 引用、无新增 `Overfull`。

---

# ApproxOperator.jl 编程规则（`program/`）

本文件同时约束本项目对 `~/code/ApproxOperator.jl` 的使用和扩展。
新增/修改 Julia 程序时，以下规则优先于通用 Julia 习惯。

## 1. 代码位置与模块

- 包路径固定为 `~/code/ApproxOperator.jl`。
- 弱形式算子写入 `src/operation/` 下对应的物理子模块，不直接写进主模块。
- 子模块写法：

```julia
module KirchhoffPlate
using ..ApproxOperator: AbstractElement
# ...
end
```

- 子模块不 `export` 算子；程序中使用：

```julia
import ApproxOperator.KirchhoffPlate: ∫κᵢⱼMᵢⱼdΩ, ∫wqdΩ
```

- 新增源文件必须在 `src/ApproxOperator.jl` 中 `include`。

## 2. 算子签名

- 单元级算子：

```julia
function ∫...dΩ(ap::T, k::AbstractMatrix{Float64};
                 D::Real, ν::Real) where T<:AbstractElement
```

- 材料常数一律用关键字参数；不得再使用旧的全局 `op` 或 `op.D` 写法。
- 边界值、精确场、载荷等由 `prescribe!` 写入积分点，算子从 `ξ` 读取，例如 `ξ.g`、`ξ.θ`、`ξ.q`。
- 双单元集耦合算子：

```julia
function ∫...dΩ(aʷ::T, aᵐ::S, k::AbstractMatrix{Float64}) where
        {T<:AbstractElement,S<:AbstractElement}
```

调用：

```julia
(∫...dΩ => (w_elements, M_elements))(k)
```

- 每个混合算子只填一个独立编号的块矩阵；转置块由主程序 `k'` 生成，不得在算子中同时接收 `k_wM` 和 `k_Mw`。
- **边界算子**：边界变量 `w` 必须放在参数最后，例如
  `∫VwdΓ(aᵐ, aʷ, k)`。
  `k` 的行是 M 自由度（`3I-2`、`3I-1`、`3I`），列是 w 自由度 `J`；
  `k` 尺寸为 `(3nM, nw)`，转置块为 `k'`。
  不得写成 `(aʷ, aᵐ, k)`，也不得写成 `k[I, 3J-2]`。
- 同一算子用重载区分两个积分区域：
  - 不带 `f` 的 `(aᵐ, aʷ, k)` 服务每个单元的 `Γ_C`，用全边界 `Γ` 的符号；
  - 带 `f` 的 `(aᵐ, aʷ, k, f)` 服务强制边界 `Γ_w`、`Γ_θ`、`c_w`，用强制边界的符号。
  两套符号相反，不得混用，也不再加 `side` 关键字或另起算子名。

## 3. 自由度约定

- 挠度 `w`：一个自由度，节点编号为 `I`。
- 弯矩 `M`：每个 piecewise 节点三个自由度：

```text
M₁₁ -> 3I-2
M₂₂ -> 3I-1
M₁₂ -> 3I
```

- 域内耦合 `∫∇²MwdΩ(aᵐ, aʷ, k)` 填 **row M, column w**：
  `Kᵐʷ`，尺寸 `(3nM, nw)`。
- 边界算子 `∫VwdΓ`、`∫MθdΓ`、`wP` 也填 **row M, column w**：
  同一个 `Kᵐʷ`。
- 弯矩柔度 `∫MC⁻¹MdΩ(elements_m, kᵐᵐ; D, ν)` 填 `kᵐᵐ`，尺寸 `(3nM, 3nM)`。
- 混合系统不 Schur 压缩，主程序组装：

```julia
K = [zeros(nw, nw)  Kᵐʷ'
     Kᵐʷ            kᵐᵐ]

f = [fʷ; fᵐ]
```

- `Kᵐʷ` 由域内耦合和边界算子共同组装，`kᵐᵐ` 由弯矩柔度算子组装。

## 4. M 场与 RKGSI

- 弯矩场 `M_{ij}` 用 `PiecewisePolynomial` 独立离散；三个分量只是不同自由度，共用同一组 piecewise 形函数。
- 域内耦合使用：

```text
∫ M_ij,ij w dΩ
```

即二阶导作用在 M 的 piecewise 形函数上，`w` 取形函数值 `Nʷ`。
- 弯矩柔度块 `∫ δM : C⁻¹ : M` 用 M 的形函数值。
- 边界 `V(M)` 用 M 的一阶导；`Γθ` 的 `Mnn`、角点 `Mns` 用 M 的形函数值。
- 角点项需要 `ξ.Δ`；默认 Gauss 积分点不含端点，`Δ` 全为零，需要含端点的边界积分规则。

## 5. 命名规则

混合块中统一使用：

```text
Nʷ, Nᵐ        形函数值
aʷ, aᵐ        单元集
ξʷ, ξᵐ        积分点
B₁₁, B₁₂, B₂₂ 二阶导
Bʷ₁, Bʷ₂      一阶导
```

不得再写：

```text
Nw, NM, aM, aw, ξM, ξw, Bxx, Bxy, Byy
```

数学算子名保留 Unicode 风格，例如 `∫κᵢⱼMᵢⱼdΩ`、`∫∇²MwdΩ`。

## 6. 形函数 setter

统一使用：

```julia
set𝝭!
set∇𝝭!
set∇²𝝭!
set∇̂³𝝭!
```

- 三阶导一律用 `set∇̂³𝝭!`。
- 每个算子的 docstring 必须写明它需要的 setter，例如：

```text
set∇²𝝭! — reads ξ[:∂²𝝭∂x²], ξ[:∂²𝝭∂x∂y], ξ[:∂²𝝭∂y²], ξ.𝑤
```

## 7. 前处理

- 边界单元需要法向时：

```julia
getElements(nodes, entities["Γᵍ"], type, order, sp;
            normal = true, delta = true)
```

- `delta = true` 生成 `ξ.Δ`；`Δ = -1`、`0`、`+1` 分别对应边界单元局部坐标的端点与内部点。
- meshfree 边界积分按 Gmsh entity 逐条边组装，避免多 entity 物理组一次性读取时的重复生成问题。

## 8. 扩展与检查规则

- 新增算子前先检查 `src/operation/` 和 `docs/src/operators.md` 中是否已有同名或等价算子。
- 找不到所需算子时先向作者报告，不得私自新增弱形式算子或改名。
- 修改后至少执行：

```bash
julia -e 'using ApproxOperator; println("loaded")'
```

- 对小单元做冒烟测试，检查块矩阵尺寸、转置关系和有限性；不能只靠“能加载”判断正确。
- 不批量重命名无关代码，不改动未涉及的物理模块。
