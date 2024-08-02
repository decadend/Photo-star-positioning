import re
import datetime
import math

# 将时角和赤纬转换为经度和纬度的函数
def calculate_lat_lon(ha_h, ha_m, ha_s, dec_d, dec_m, dec_s):
    # 计算纬度
    lat = math.radians(dec_d + dec_m / 60 + dec_s / 3600)

    # 计算经度
    lon = 2 * math.pi - math.radians((ha_h + ha_m / 60 + ha_s / 3600) * 15)

    return lat, lon


# 计算法向量的函数
def calculate_normal_vector(lat, lon):
    a = math.cos(lon) * math.cos(lat)
    b = math.sin(lon) * math.cos(lat)
    c = math.sin(lat)
    return a, b, c


# 计算两个法向量之间的夹角余弦值
def calculate_cosine_of_angle(vector1, vector2):
    dot_product = sum(a * b for a, b in zip(vector1, vector2))
    return dot_product


# 计算实际夹角和焦距的函数
def calculate_actual_angle_and_focal_length(coords1, coords2, cosine_value):
    x1, y1 = coords1
    x2, y2 = coords2

    # 计算 d, e, f 和 g
    d = cosine_value ** 2 - 1
    e = ((x1 ** 2 + y1 ** 2 + x2 ** 2 + y2 ** 2) * cosine_value ** 2 - 2 * (x1 * x2 + y1 * y2))
    f = ((x1 ** 2 + y1 ** 2) * (x2 ** 2 + y2 ** 2) * cosine_value ** 2 - (x1 * x2 + y1 * y2) ** 2)
    g = e ** 2 - 4 * d * f

    # 输出每一步的值
    print(f"d: {d}")
    print(f"e: {e}")
    print(f"f: {f}")
    print(f"g: {g}")

    # 计算焦距
    if d != 0 and g >= 0:
        focal_length = math.sqrt((-math.sqrt(g) - e) / (2 * d))
    else:
        focal_length = float('nan')  # 处理 d 为 0 或 g 小于 0 的情况

    return g, focal_length


# 定义时角和赤纬
data = [
    (13, 48, 53.56, +54, 47, 49.2),
    (13, 25, 18.01, +49, 11, 25.9),
    (14, 16, 37.74, +38, 11, 11.9),
    (15, 26, 27.59, +47, 38, 36.0),
    (16, 2, 45.68, +44, 21, 55.6)
]
# 定义天顶坐标
zenith_x, zenith_y = 80.9139979103, -1166.4665302339438
# 定义平面坐标
coordinates = [
    (-215.5111361369099, -36.4210622433514),
    (-111.0612502625295, 35.6892987765508),
    (94.2976296443091, -139.6687807287531),
    (-132.8747956681829, -354.5500668256116),
    (-117.9846485353764, -511.469074473517)
]

# 计算所有法向量
vectors = []
for ha_h, ha_m, ha_s, dec_d, dec_m, dec_s in data:
    lat, lon = calculate_lat_lon(ha_h, ha_m, ha_s, dec_d, dec_m, dec_s)
    vector = calculate_normal_vector(lat, lon)
    vectors.append(vector)
# 输出每个法向量的值
for i, vector in enumerate(vectors):
    print(f"法向量 {i + 1}: {vector}")
# 存储所有焦距和仰角的 sin 值
focal_lengths = []
sine_angles = []  # 存储 sin(仰角) 的值

# 计算并输出每对法向量之间的夹角余弦值、实际夹角和焦距
for i in range(len(vectors)):
    for j in range(i + 1, len(vectors)):
        cosine_of_angle = calculate_cosine_of_angle(vectors[i], vectors[j])
        print(f"法向量 {i + 1} 和 {j + 1} 之间的夹角余弦值: {cosine_of_angle}")
        g, focal_length = calculate_actual_angle_and_focal_length(coordinates[i], coordinates[j], cosine_of_angle)
        if not math.isnan(focal_length):  # 忽略无效的焦距
            focal_lengths.append(focal_length)
        print(f"法向量 {i + 1} 和 {j + 1} 之间的实际夹角计算值: {g}")
        print(f"法向量 {i + 1} 和 {j + 1} 之间的焦距: {focal_length}\n")

# 计算平均值
if focal_lengths:
    normal_average = sum(focal_lengths) / len(focal_lengths)

print(f"焦距的平均值: {normal_average}")



# 计算天顶坐标的膜长
zenith_length = math.sqrt(zenith_x ** 2 + zenith_y ** 2 + normal_average ** 2)
print(f"天顶坐标的膜长: {zenith_length}")

# 计算每颗星星的膜长
star_lengths = []
for (x, y) in coordinates:
    star_length = math.sqrt(x ** 2 + y ** 2 + normal_average ** 2)
    star_lengths.append(star_length)
    print(f"星星坐标 ({x}, {y}) 的膜长: {star_length}")

# 计算每颗星星的仰角（以弧度为单位）并计算 sin(仰角)
angles = []  # 存储每颗星星的仰角弧度
sine_angles = []  # 存储 sin(仰角) 的值
for i, star_length in enumerate(star_lengths):
    # 使用公式计算仰角
    numerator = zenith_x * coordinates[i][0] + zenith_y * coordinates[i][1] + normal_average ** 2
    denominator = zenith_length * star_length
    angle = math.pi / 2 - math.acos(numerator / denominator)
    angles.append(angle)
    sin_angle = math.sin(angle)
    sine_angles.append(sin_angle)
    print(f"第 {i + 1} 颗星星的仰角 (弧度): {angle}")
    print(f"第 {i + 1} 颗星星的 sin(仰角): {sin_angle}")

# 打开文件以写入经纬度
with open('/Users/Zhuanz/Downloads/坐标纬度加二组4.txt', 'w') as file:
    # 计算每对法向量之间的公式值
    print("\n计算每对法向量之间的公式值:")
    for i in range(len(vectors)):
        for j in range(i + 1, len(vectors)):
            vector_a = vectors[i]
            vector_b = vectors[j]

            # 计算公式1、公式2和公式3
            cross_product_x = vector_b[1] * vector_a[2] - vector_b[2] * vector_a[1]
            cross_product_y = vector_b[2] * vector_a[0] - vector_b[0] * vector_a[2]
            cross_product_z = vector_b[0] * vector_a[1] - vector_b[1] * vector_a[0]
            formula1_value = cross_product_x ** 2 + cross_product_y ** 2 + cross_product_z ** 2
            print(f"法向量 {i + 1} 和 {j + 1} 之间的公式1的计算值: {formula1_value}")

            # 计算公式2的计算步骤
            sin_angle_i = sine_angles[i]
            sin_angle_j = sine_angles[j]

            terma = vector_b[0] * vector_a[2] - vector_a[0] * vector_b[2]
            termb = vector_a[2] * sin_angle_j - vector_b[2] * sin_angle_i
            termc = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0]
            termd = vector_b[1] * sin_angle_i - vector_a[1] * sin_angle_j

            term1 = 2 * (terma * termb + termc * termd)
            print(f"法向量 {i + 1} 和 {j + 1} 之间的公式2的计算值: {term1}")

            # 计算公式3的计算步骤
            qiugenga = vector_a[2] * sin_angle_j - vector_b[2] * sin_angle_i
            qiugengb = vector_b[1] * sin_angle_i - vector_a[1] * sin_angle_j
            qiugengc = vector_a[1] * vector_b[2] - vector_b[1] * vector_a[2]
            qiugengd = qiugenga**2 + qiugengb**2 - qiugengc**2

            print(f"法向量 {i + 1} 和 {j + 1} 之间的公式3的计算值: {qiugengd}")

            # 计算 bb^2 - 4 * aa * cc
            aa = formula1_value
            bb = term1
            cc = qiugengd
            dd = bb**2 - 4 * aa * cc

            # 输出结果到控制台
            print(f"法向量 {i + 1} 和 {j + 1} 之间的 bb^2 - 4 * aa * cc 的计算值: {dd}")

            # 计算一元二次方程的根
            if aa != 0 and dd >= 0:  # 确保 a 不为零且判别式非负
                sqrt_dd = math.sqrt(dd)

                # 计算 x1 和 x2
                x1 = (-bb + sqrt_dd) / (2 * aa)
                x2 = (-bb - sqrt_dd) / (2 * aa)
                print(f"法向量 {i + 1} 和 {j + 1} 之间的 x1 的计算值: {x1}")
                print(f"法向量 {i + 1} 和 {j + 1} 之间的 x2 的计算值: {x2}")

                # 计算 y1 和 y2 的值
                terme = vector_b[0] * vector_a[2] - vector_a[0] * vector_b[2]
                termb = vector_a[2] * sin_angle_j - vector_b[2] * sin_angle_i
                termf = vector_a[1] * vector_b[2] - vector_b[1] * vector_a[2]

                if termf != 0:
                    # 计算 y1 对应 x1 的值
                    y1_x1 = (terme * x1 + termb) / termf
                    # 计算 y2 对应 x2 的值
                    y2_x2 = (terme * x2 + termb) / termf
                    print(f"法向量 {i + 1} 和 {j + 1} 之间的 y1 对应 x1 的计算值: {y1_x1}")
                    print(f"法向量 {i + 1} 和 {j + 1} 之间的 y2 对应 x2 的计算值: {y2_x2}")

                    # 计算 z1 和 z2 的值
                    termg = vector_a[0] * vector_b[1] - vector_b[0] * vector_a[1]
                    termh = vector_b[1] * sin_angle_i - vector_a[1] * sin_angle_j
                    termini = vector_b[1] * vector_a[2] - vector_a[1] * vector_b[2]

                    if termini != 0:
                        # 计算 z1 对应 x1 的值
                        z1_x1 = (termg * x1 + termh) / termini
                        # 计算 z2 对应 x2 的值
                        z2_x2 = (termg * x2 + termh) / termini
                        print(f"法向量 {i + 1} 和 {j + 1} 之间的 z1 对应 x1 的计算值: {z1_x1}")
                        print(f"法向量 {i + 1} 和 {j + 1} 之间的 z2 对应 x2 的计算值: {z2_x2}")

                        # 转换为经纬度
                        def convert_to_lat_lon(x, y, z):
                            longitude = math.degrees(math.atan2(-y, -x))  # 使用 ATAN2(-y, -x) 计算经度
                            latitude = math.degrees(math.asin(z / math.sqrt(x**2 + y**2 + z**2)))  # 使用 ASIN(z) 计算纬度
                            return latitude, longitude

                        lat1_x1, lon1_x1 = convert_to_lat_lon(x1, y1_x1, z1_x1)
                        lat2_x2, lon2_x2 = convert_to_lat_lon(x2, y2_x2, z2_x2)

                        # 将经纬度以纬度，经度的格式写入文件中
                        file.write(f"{lat1_x1:.6f},{lon1_x1:.6f}\n")
                        file.write(f"{lat2_x2:.6f},{lon2_x2:.6f}\n")
                    else:
                        file.write(f"法向量 {i + 1} 和 {j + 1} 之间的 termini 为零，不能计算 z1 和 z2。\n")
                else:
                    file.write(f"法向量 {i + 1} 和 {j + 1} 之间的 termf 为零，不能计算 y1 和 y2。\n")
            elif dd < 0:
                file.write(f"法向量 {i + 1} 和 {j + 1} 之间的判别式小于0，无实数根，跳过 x1, x2, y1, y2 和 z1, z2 计算\n")
            else:
                file.write(f"法向量 {i + 1} 和 {j + 1} 之间的 a 不能为零，请检查输入值。\n")

print("经纬度已导出到 '/Users/Zhuanz/Downloads/坐标.txt' 文件中。")