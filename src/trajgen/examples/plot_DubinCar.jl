using Plots
using Measures

"""
    plot_dubins_trajectory_final(...)

根据指定的坐标系定义，绘制 Dubins Car 的轨迹。

# 坐标系定义
- x: 横向坐标
- y: 纵向坐标
- theta: 车辆纵向指向与 Y 轴正向的夹角。指向 X 轴正向为正角度（即从Y轴开始逆时针为正）。

# 参数
- `trajectory`: 包含 [x, y, theta] 状态向量的 Vector{Vector{Float64}}。
- `N_lanes`: 车道数量 (默认为 3)。
- `lane_width`: 每条车道的宽度 (默认为 3.75 米)。
- `road_length`: 道路的纵向长度 (默认为 100.0 米)。
- `vehicle_length`: 模拟车辆的长度 (默认为 3.0m)。
- `vehicle_width`: 模拟车辆的宽度 (默认为 1.5m)。
- `orientation_line_length`: 朝向指示线的长度 (默认为 2.5m)。
"""
function plot_AutoTrjPbm_DubinCar(trajectory::Vector{Vector{Float64}};
                                      N_lanes::Int=3,
                                      lane_width::Float64=3.75,
                                      road_length::Float64=100.0,
                                      vehicle_length::Float64=3.0,
                                      vehicle_width::Float64=1.5,
                                      orientation_line_length::Float64=2.5,
                                      arrow_frequency::Int=1)
    
    # 创建一个 Plot 对象
    p = plot(aspect_ratio=:equal,
             xlabel="X position (m) - Lateral",
             ylabel="Y position (m) - Longitudinal",
             title="Dubins Car Trajectory",
             legend=:none, # <-- 修改: 图例放置于右下角
             size=(600, 750),     # <-- 修改: 按比例缩小图像尺寸
             margin=10mm)

    # 1. 绘制车道
    y_min_road, y_max_road = 0.0, road_length * 1.05

    for i = 0:N_lanes
        plot!(p, [i * lane_width, i * lane_width], [y_min_road, y_max_road],
              line=:dash, color=:gray, label=(i == 0 ? "Lane Boundary" : nothing))
    end
    
    for i = 1:N_lanes
        plot!(p, [(i - 0.5) * lane_width, (i - 0.5) * lane_width], [y_min_road, y_max_road],
              line=:dot, color=:lightgray, label=(i == 1 ? "Lane Center" : nothing))
    end

    # 2. 绘制车辆主体矩形和朝向线段
    for i in 1:arrow_frequency:length(trajectory)
        x, y, theta = trajectory[i]
        
        phi = pi/2 - theta

        # 车辆主体矩形
        half_length = vehicle_length / 2
        half_width = vehicle_width / 2
        
        corners_local = [
             half_length  half_width;
             half_length -half_width;
            -half_length -half_width;
            -half_length  half_width;
             half_length  half_width;
        ]

        R = [cos(phi) -sin(phi); sin(phi) cos(phi)]
        
        corners_global_x = zeros(size(corners_local, 1))
        corners_global_y = zeros(size(corners_local, 1))
        for j = 1:size(corners_local, 1)
            rotated_point = R * corners_local[j, :]
            corners_global_x[j] = x + rotated_point[1]
            corners_global_y[j] = y + rotated_point[2]
        end
        
        plot!(p, corners_global_x, corners_global_y,
              line=(1.5, :black, :solid),
              fillcolor=:gray, fillalpha=0.8,
              label=(i == 1 ? "Vehicle Body" : nothing))

        # 车辆朝向指示线
        dx_orient = orientation_line_length * cos(phi)
        dy_orient = orientation_line_length * sin(phi)
        
        plot!(p, [x, x + dx_orient], [y, y + dy_orient],
              line=(2, :red, :solid),
              label=(i == 1 ? "Orientation" : nothing))
    end

    # 3. 标记起点和终点
    start_state = trajectory[1]
    end_state = trajectory[end]

    scatter!(p, [start_state[1]], [start_state[2]],
             marker=(:circle, 8, :green, :fill),
             label="Start: ($(round(start_state[1], digits=2)), $(round(start_state[2], digits=2)))")

    scatter!(p, [end_state[1]], [end_state[2]],
             marker=(:star5, 10, :purple, :fill),
             label="End: ($(round(end_state[1], digits=2)), $(round(end_state[2], digits=2)))")

    # 4. 设置固定的坐标轴范围
    xlims!(p, -3.75, 4 * 3.75)
    ylims!(p, -5, road_length + 5)
             
    return p
end
