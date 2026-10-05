# GNSS观测分析后端212

当前仓库：https://github.com/huangjie666777-ux/gnss-reflectometry-212

纯后端服务：上传一份未压缩 RINEX 3.04 观测文件和一份 SP3-c 精密星历：

- /position — 逐历元独立解算接收机 ECEF 位置与钟差（C1C 单点定位）；
- /tec — GPS 双频（C1C/C2W/L1C/L2W）电离层 TEC 监测，输出逐星逐历元
  斜向/垂直 TEC、弧编号、穿刺点地心经纬度与排除原因。
- /reflect — 单频 GNSS 反射测高（GNSS-IR）：仅用 S1C 信噪比弧段反演
  反射面高度，输出水位（Z - H）、幅值、残差 RMS 与高度扫描曲线。

## 环境

- Python 3.10.12 / FastAPI 0.115.12 / NumPy 2.2.6（全部在 .venv 中）

## 运行

    .venv/bin/python -m uvicorn main:app --host 127.0.0.1 --port 8000

## 使用

    curl -s -F "rinex=@examples/obs.rnx" -F "sp3=@examples/eph.sp3" \
         http://127.0.0.1:8000/position | python3 -m json.tool

TEC 监测（站点 WGS84 经纬高 + 逐星合并码偏差 JSON，单位 ns）：

    curl -s -F "rinex=@examples/obs.rnx" -F "sp3=@examples/eph.sp3" \
         -F "station_lat_deg=30.0" -F "station_lon_deg=114.0" \
         -F "station_height_m=50.0" \
         -F "biases=$(cat examples/biases.json)" \
         http://127.0.0.1:8000/tec | python3 -m json.tool

反射测高（站点 WGS84 经纬高、天线相位中心相对水尺零点高程 Z，
以及正递增的反射高度搜索上下界与正步长，单位米）：

    curl -s -F "rinex=@examples/obs_reflect.rnx" -F "sp3=@examples/eph_reflect.sp3" \
         -F "station_lat_deg=30.0" -F "station_lon_deg=114.0" \
         -F "station_height_m=50.0" -F "antenna_height_z_m=12.0" \
         -F "reflector_min_m=1.0" -F "reflector_max_m=15.0" \
         -F "reflector_step_m=0.01" \
         http://127.0.0.1:8000/reflect | python3 -m json.tool

返回每个历元：状态与失败原因、ECEF 米坐标、WGS84 经纬度与椭球高、
接收机钟差（秒）、使用/排除的卫星（含排除原因）、逐星残差（米）与 RMS。

## 生成合成示例数据

    .venv/bin/python examples/make_synthetic.py   # 写 examples/obs.rnx, examples/eph.sp3
    .venv/bin/python examples/make_reflect_sample.py  # 写 examples/obs_reflect.rnx, examples/eph_reflect.sp3

反射测高样例含已知水位：Z = 12.0 m，真实反射高度 H = 8.0 m，
即已知水位 Z - H = 4.0 m（G01 升弧 / G02 落弧 / G03 含 S1C 缺口断弧 /
G04 升落反转分弧 / G05 全空 S1C 记录）。观测文件头部仅声明 S1C，
用于验证"允许观测头只有 S1C"。

## 测试

    .venv/bin/python -m pytest tests -q

## 编译检查

    .venv/bin/python -m compileall -q gnss_reflect212 examples main.py

## 处理规则与限制

- 仅 GPS 时间系统、GPS 卫星、C1C 伪距、正常历元（flag 0）；
  接收机钟改正预应用（RCV CLOCK OFFS APPL != 0）、事件历元、非 3.04 版本、
  压缩文件等明确拒绝并指明文件与行号。
- 按头部观测类型顺序解析固定宽度字段，兼容 SYS / # / OBS TYPES 续行；
  空白或零伪距视为缺测。
- 校验日期、卫星身份、历元数量（上限 100）、非有限值与截断。
- SP3：位置 km→m，钟差 µs→s；零坐标与 999999.999999 缺失钟差被识别。
- 卫星位置用连续 8 节点 Lagrange 插值，钟差用相邻节点线性插值；
  不跨缺测节点、不外推，不可用的卫星给出排除原因。
- 每历元独立迭代最小二乘（位置 + 钟差），由伪距与卫星钟差推发射时刻，
  补偿信号飞行期间地球自转（Sagnac）；至少 4 颗有效卫星且几何满秩才求解；
  头部近似坐标只作初值，不沿用上一历元结果。
- 不足 4 星、秩亏或不收敛的历元返回原因，其余历元继续。

## 反射测高解算规则

- 模型假设：反射面为静止水平平面，不计大气折射；多路径振荡相位为
  4*pi*H*sin(仰角)/lambda_L1，H 为天线相位中心相对反射面的垂直高度，
  水位 = Z - H（Z 为天线相位中心相对水尺零点高程）。
- 仅 GPS 时间系统、正常历元（flag 0）与 S1C 观测量，最多 100 历元；
  观测头允许只有 S1C（定位与 TEC 端点各自检查所需观测类型并明确拒绝）。
- 卫星位置复用定位链路的连续 8 节点 Lagrange 插值（接收时刻取值），
  不跨缺失节点、不外推；仰角由站点局部 ENU 坐标计算，仅取 5-25 度。
- 分弧：逐星按时间排列，缺 S1C、无星历（含插值跨缺失节点）、
  相邻间隔 >120 s、仰角超出 5-25 度过滤窗及升落方向反转均断弧，
  弧不跨过滤缺口；被排除的历元保留在记录中并注明原因。
- 弧质量门限：至少 12 个样本且仰角跨度至少 5 度，否则标 failed。
- 反演：S1C 按 10^(S1C/20) 转线性幅度，在 sin(仰角) 上拟合并减去
  二次趋势；因 sin(仰角) 非均匀采样不做 FFT，改为在规则网格
  （上界 5001 格）扫描 H，对每个 H 最小二乘拟合
  cos/sin(4*pi*H*sin(仰角)/lambda_L1) 加常数项，取残差平方和最小者，
  同值取较小 H；秩亏或去趋势后无剩余波动的弧标 failed，
  最优落在网格边界的弧给出警告。
- 每弧返回：起止时间、升落方向、样本来源（S1C）、反射高度 H、
  水位 Z - H、幅值、残差 RMS、警告及高度扫描曲线（H 与对应 RSS）。

## 修复记录

- 全空卫星记录不再消失：观测行所有字段为空时仍保留该卫星记录，
  下游（TEC/反射）逐历元给出排除原因。
- 星历插值不再跨缺失节点：8 节点窗口与钟差相邻节点均校验时间
  间隔连续性，跨缺口时返回排除原因而非静默插值。

## TEC 解算规则

- 观测组合：码几何无关量 GF_P = (P2 - P1) - c·DCB，载波几何无关量
  GF_L = L1·λ1 - L2·λ2（周转米）；两者均等于 40.3·STEC·(1/f2² - 1/f1²)，
  载波多一个弧内常数模糊度偏移。
- 偏差符号约定：biases 给出逐星合并 P1-P2 码偏差（DCB，ns），即码观测量
  P2 相对 P1 多含 +c·DCB 硬件延迟，因此从 (P2 - P1) 中**减去** c·DCB；
  缺偏差的卫星整星标 invalid，不静默丢弃。
- 分弧：四量缺一、L1C/L2W 的 LLI 非零、相邻历元间隔 >120 s 或载波 GF
  跳变 >1 m 均断弧；各弧独立定级（不共享偏移），以弧内
  mean(GF_P - GF_L) 为偏移加到载波上；不足 3 个完整样本的弧标 invalid。
  缺测历元在输出中保留并注明原因，不静默消失。
- 斜向 TEC = GF_L(定级后) / [40.3·(1/f2² - 1/f1²)]，单位 TECU，保留负值。
- 薄壳近似：站星射线与半径 6821 km 地心球壳求交得穿刺点（输出地心
  经纬度）；垂直 TEC = 斜向 TEC × cos(射线与壳面法向夹角)。
  站点径向取为天顶，仰角 <10° 的历元标 excluded 并注明原因。
- 卫星位置复用定位链路的 8 节点 Lagrange 插值（接收历元处取值），
  不跨缺失节点、不外推；插值不可用的历元标 excluded。
- 站点坐标与偏差做有限值/范围校验；未支持的修正
  （如 RCV CLOCK OFFS APPL ≠ 0）明确拒绝。

## 精度边界

仅 C1C 伪距单点定位，无电离层/对流层/载波/多路径改正，
合成数据下亚米级；实际观测精度为米级，受大气延迟与观测噪声主导。

## 模块划分

- gnss_reflect212/rinex.py — RINEX 3.04 解析与校验
- gnss_reflect212/sp3.py — SP3-c 解析与单位换算
- gnss_reflect212/interp.py — 位置 Lagrange / 钟差线性插值
- gnss_reflect212/solver.py — 逐历元最小二乘定位
- gnss_reflect212/tec.py — 双频几何无关组合、分弧定级、薄壳映射与穿刺点
- gnss_reflect212/reflect.py — S1C 弧划分、去趋势与网格扫描反射高度反演
- gnss_reflect212/geodesy.py — WGS84 坐标转换与常数
- gnss_reflect212/main.py — FastAPI 入口
