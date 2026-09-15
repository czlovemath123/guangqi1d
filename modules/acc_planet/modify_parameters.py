import sys

p = sys.argv[1:]
if len(p) < 4:
    print("Parameter not enough!")
    sys.exit(1)
elif len(p) > 4:
    print("Too many parameters!")
    sys.exit(1)
else:
    mdot = p[0]
    mp = p[1]
    rp = p[2]
    t_kh = p[3] 

    #处理problem.data
    # 读取文件内容
    with open('problem.data', 'r') as file:
        data = file.readlines()

    # 遍历文件中的每一行并查找需要替换的内容
    for i in range(len(data)):
        if 'acc_rate=' in data[i]:
            data[i] = f"    acc_rate={mdot}d-3               !accretion rate in earth mass per year\n"
        elif 'm_planet=' in data[i]:  # 查找
            data[i] = f"    m_planet={mp}d0                  !planetary mass in Jupiter mass\n"
        elif 'kh_timescale=' in data[i]:
            data[i] = f"    kh_timescale={t_kh}e5             !in years\n"

    # 将修改后的内容写回文件
    with open('problem.data', 'w') as file:
        file.writelines(data)
    print("problem.data modify successfully")
    
    #处理global.data
    # 读取文件内容
    with open('global.data', 'r') as file:
        gdata = file.readlines()
    
    # 遍历文件中的每一行并查找需要替换的内容
    for j in range(len(gdata)):
        if 'lengthscale = ' in gdata[j]:
            gdata[j] = f"    lengthscale = {rp}d9\n"

    # 将修改后的内容写回文件
    with open('global.data', 'w') as file:
        file.writelines(gdata)
    print("global.data modify successfully")
    sys.exit(0)

