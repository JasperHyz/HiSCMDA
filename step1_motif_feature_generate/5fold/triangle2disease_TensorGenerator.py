# coding=utf-8
import sys
import getopt

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

input_file = 'triangle2disease_motif_SEED_{0}_FOLD_{1}'.format(SEED, FOLD)
output_file = 'triangle2disease_step1_SEED_{0}_FOLD_{1}'.format(SEED, FOLD)
dic_file_template = '{0}_dic'


f = open(input_file)
tensor = [line.strip('\n').split("\t") for line in f]
id_dic = {}
next_id = 1

for i in range(0, len(tensor)):
    for j in range(0, len(tensor[i])):
        if id_dic.get(tensor[i][j]) is None:
            id_dic[tensor[i][j]] = next_id
            next_id = next_id + 1

content = ''
for key, value in id_dic.items():
    content += ("" if len(content) == 0 else "\n") + "{0} {1}".format(key, value)

dic_file_path = dic_file_template.format(output_file)
dic_file = open(dic_file_path, 'w')
dic_file.write(content)
dic_file.close()

content = ''
for i in range(0, len(tensor)):
    # tensor[i].sort()
    mapped_tensor = []
    for j in range(0, len(tensor[i])):
        mapped_tensor.append(id_dic.get(tensor[i][j]))
    max_index = mapped_tensor.index(max(mapped_tensor))

    new_tensor = []
    for j in range(max_index, len(mapped_tensor)):
        new_tensor.append(mapped_tensor[j])

    for j in range(0, max_index):
        new_tensor.append(mapped_tensor[j])

    content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(new_tensor[0],
                                                                              new_tensor[1],
                                                                              new_tensor[2], 1)
    content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(new_tensor[0],
                                                                              new_tensor[2],
                                                                              new_tensor[1], 1)
    content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(new_tensor[1],
                                                                              new_tensor[0],
                                                                              new_tensor[2], 1)
    content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(new_tensor[1],
                                                                              new_tensor[2],
                                                                              new_tensor[0], 1)
    content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(new_tensor[2],
                                                                              new_tensor[1],
                                                                              new_tensor[0], 1)
    content += ("" if len(content) == 0 else "\n") + '{0} {1} {2} {3}'.format(new_tensor[2],
                                                                              new_tensor[0],
                                                                              new_tensor[1], 1)

new_tensor_file = open(output_file, 'w')
new_tensor_file.write(content)
new_tensor_file.close()

