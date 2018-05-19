unstaged_results_file = "E:\\MSC\project\\Deliverables\\results-unstagged.txt"
staged_results_file = "E:\\MSC\\project\\Deliverables\\results-stagged.txt"

unstaged_content_array = []
staged_content_array = []

output_file = "E:\\MSC\\project\\Deliverables\\comparison.txt"

unstaged_contents = []
staged_contents = []


with open(unstaged_results_file) as f:
    unstaged_contents = f.readlines()

with open(staged_results_file) as f:
    staged_contents = f.readlines()

for i in range(len(unstaged_contents)):
    content = unstaged_contents[i]
    content = content.rstrip()
    content = content.split()
    num_content = []
    for j in content:
        if "nan" in j:
            num_content.append(0)
        else:
            num_content.append(float(j))
    unstaged_content_array.append(num_content)



for i in range(len(staged_contents)):
    content = staged_contents[i]
    content = content.rstrip()
    content = content.split()
    num_content = []
    for j in content:
        if "nan" in j:
            num_content.append(0)
        else:
            num_content.append(float(j))
    staged_content_array.append(num_content)



f = open(output_file, "w")
for i in range(len(unstaged_content_array)):
    for j in range(len(unstaged_content_array[i])):
        if(unstaged_content_array[i][j]!=staged_content_array[i][j]):
            f.write("Mismatch in "+str(i)+","+str(j)+"\n")
            f.write("Difference is "+str(unstaged_content_array[i][j]-staged_content_array[i][j])+str("\n"))
            f.write("Unstagged-Value:"+str(unstaged_content_array[i][j])+str("\n"))
            f.write("Stagged-Value:" + str(staged_content_array[i][j]) + str("\n"))

