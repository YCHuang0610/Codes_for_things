import os
import shutil
#用来把所有子文件夹中的文件挪出来的小程序，用来汇总照片图片。
#把该程序放在目录中双击就行
#递归寻找文件
def findfile(working_path, file_list):
    for file_or_dir in os.listdir(working_path):
        file_or_dir_path = os.path.join(working_path, file_or_dir)
        if os.path.isfile(file_or_dir_path):
            file_list.append(file_or_dir_path)
        else:
            findfile(file_or_dir_path, file_list)
            
#移动文件列表中所有文件到目标目录并重命名
def move(des_dir, file_list):
    for num, file in enumerate(file_list, start=1):
        file_type = os.path.splitext(file)[1]
        shutil.move(file, str(num)+file_type)

#去除指定文件名
def dele(file_list, name):
    i = 0
    while i < len(file_list):
        file_name = os.path.basename(file_list[i])
        if file_name == name or file_name == os.path.basename(__file__): #除去自己
            del file_list[i] 
            i -= 1  #注意删了之后列表变短
        i += 1


if __name__=='__main__':
    working_dir = os.getcwd()
    print(working_dir)
    file_list = []
    findfile(working_dir, file_list)
    dele(file_list, ".DS_Store")
    move(working_dir, file_list)
    print("OKKK")
