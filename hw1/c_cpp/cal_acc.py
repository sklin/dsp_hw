#!/usr/bin/env python2

import sys

if len(sys.argv) < 3:
    print "[Usage] cal_acc <result_file> <answer_file>"
    sys.exit(0)

resultFile = sys.argv[1]
answerFile = sys.argv[2]

result = open(resultFile)
answer = open(answerFile)

result = [line.strip().split()[0] for line in result.readlines()]
answer = [line.strip() for line in answer.readlines()]

successCount = 0

for i in range(len(result)):
    if(result[i] == answer[i]):
        successCount += 1

print successCount/float(len(result))
