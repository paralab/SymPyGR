fw = open('temp_cuda_bssneqs.h', 'w+');
fwCpp = open('temp_bssneqs.cpp', 'w+');

roundingPower = 20;
offsetNumber = str(10 ** roundingPower);
with open('test_cuda_bssneqs.h') as f:
    lines = f.readlines()
    for line in lines:
        if not "/" in line:
            l = line.split("=");
            tempLine = "";
            if len(l) == 2:
                endIndex = l[1].find(";")
                l[1] = l[1][:endIndex]
                temp = "roundf((" + l[1] + ")*" + offsetNumber + ")/" + offsetNumber + ";"
                tempLine = l[0] + " = " + temp
            else:
                tempLine = l[0];
            fw.write("%s\n" % tempLine)
        else:
            fw.write(line)

with open('test_bssneqs.cpp') as fCpp:
    lines = fCpp.readlines()
    for line in lines:
        if not "/" in line:
            l = line.split("=");
            tempLine = "";
            if len(l) == 2:
                endIndex = l[1].find(";")
                l[1] = l[1][:endIndex]
                temp = "roundf((" + l[1] + ")*" + offsetNumber + ")/" + offsetNumber + ";"
                tempLine = l[0] + " = " + temp
            else:
                tempLine = l[0];
            fwCpp.write("%s\n" % tempLine)
        else:
            fw.write(line)
