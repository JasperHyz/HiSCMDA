# coding=utf-8
import sys
import getopt

def tensor_max_index_search_modify(tensor, _dic):
    max_value = len(_dic)
    min_value = 1
    
    for i in range(len(tensor)):
        if tensor[i][0] == max_value:
            return True
    max_key = list(_dic.items())[len(_dic)-1][0]
    min_key = list(_dic.items())[0][0]
    temp = _dic[max_key]
    _dic[max_key] = _dic[min_key]
    _dic[min_key] = temp
    return False


SEED = 0
FOLD = 0
try:
    opts, args = getopt.getopt(sys.argv[1:], "s:f:", ['SEED=', 'FOLD='])
except getopt.GetoptError:
    print('exception')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-s', '--SEED'):
        SEED = int(arg)
    if opt in ('-f', '--FOLD'):
        FOLD = int(arg)

input_file = 'mdm_motif_SEED_{0}_FOLD_{1}'.format(SEED, FOLD)
output_file = 'mdm_step1_SEED_{0}_FOLD_{1}'.format(SEED, FOLD)
dic_file_template = '{0}_dic'

f = open(input_file)
tensor = [line.strip('\n').split("\t") for line in f]
f.close()

id_dic = {}
next_id = 1
for i in range(0, len(tensor)):
    for j in range(0, len(tensor[i])):
        if id_dic.get(tensor[i][j]) is None:
            id_dic[tensor[i][j]] = next_id
            next_id = next_id + 1

tensor_temp = []
for i in range(0, len(tensor)):
    mapped_tensor = []
    for j in range(0, len(tensor[i])):
        mapped_tensor.append(id_dic.get(tensor[i][j]))
    tensor_temp.append(mapped_tensor)   


if tensor_max_index_search_modify(tensor_temp, id_dic):
    content = ''
    for key, value in id_dic.items():
        content += ("" if len(content) == 0 else "\n") + "{0} {1}".format(key, value)
    dic_file_path = dic_file_template.format(output_file)
    dic_file = open(dic_file_path, 'w')
    dic_file.write(content)
    dic_file.close()
    
    content = ''
    for i in range(0, len(tensor_temp)):
        # mmd
        content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(tensor_temp[i][0],
                                                                                  tensor_temp[i][2],
                                                                                  tensor_temp[i][1], 1)

        content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(tensor_temp[i][2],
                                                                                  tensor_temp[i][0],
                                                                                  tensor_temp[i][1], 1)
    new_tensor_file = open(output_file, 'w')
    new_tensor_file.write(content)
    new_tensor_file.close()
else:
    content = ''
    for key, value in id_dic.items():
        content += ("" if len(content) == 0 else "\n") + "{0} {1}".format(key, value)
    dic_file_path = dic_file_template.format(output_file)
    dic_file = open(dic_file_path, 'w')
    dic_file.write(content)
    dic_file.close()    
    f = open(input_file)
    tensor_origin = [line.strip('\n').split("\t") for line in f]
    f.close()
    content = ''
    for i in range(0, len(tensor_origin)):
        mapped_tensor = []
        for j in range(0, len(tensor_origin[i])):
            mapped_tensor.append(id_dic.get(tensor_origin[i][j]))

        content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(mapped_tensor[0],
                                                                                  mapped_tensor[2],
                                                                                  mapped_tensor[1], 1)

        content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(mapped_tensor[2],
                                                                                  mapped_tensor[0],
                                                                                  mapped_tensor[1], 1)
    new_tensor_file = open(output_file, 'w')
    new_tensor_file.write(content)
    new_tensor_file.close()
