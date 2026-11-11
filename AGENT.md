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

## 编译（唯一事实来源：Helix 的 texlab 配置）

编译方式由 `~/.config/helix/languages.toml` 中 `[language-server.texlab.config.texlab.build]` 定义，等价命令为：

```bash
cd manuscript && latexmk -xelatex -auxdir=.aux -interaction=nonstopmode -file-line-error -synctex=1 main.tex
```

- 引擎固定 XeLaTeX，参数固定为上面这一组；不得再加入 `-pdf`、`-draftmode`、`-outdir`、`--no-focus` 等开关。其中 `-pdf` 会与 `-xelatex` 冲突、`-outdir` 会让 texlab 找不到文件、`--no-focus` 是 sioyek 的参数（只应出现在 `forwardSearch.args` 里，写作 `--nofocus`）。
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
