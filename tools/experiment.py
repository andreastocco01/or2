#!/usr/bin/env python3

from subprocess import run, PIPE
import itertools

# Fixed parameters
output_file = "out.txt"
executable_name = "./main"
pfflag = "--parsefriendly"
timelimit = "30"

def run_experiment(parameters) -> str:
    split_params = list(map(lambda x: x.split(), parameters))

    chained = []
    for elem in split_params:
        for elemm in elem:
            chained.append(elemm)
    command = [executable_name, pfflag, "-t", timelimit, *chained]
    print("Running {}".format(command))
    result = run(command, stdout=PIPE)
    return result.stdout.decode('utf-8')

def format_argument(key, value):
    print("key: {}, value: {}".format(key, value))
    if key != "":
        return "--{} {}".format(key, value)
    return "{}".format(value)

def generate_parameters(parameters):
    full = []

    for parameter in parameters:
        values = []
        for value in parameter[1]:
            values.append(format_argument(parameter[0], value))
        full.append(values)

    return itertools.product(*full)

def parse_result(result: str) -> (float, float):
    result = result.replace("\n ", "")
    parts = result.split(";")
    time = float(parts[0])
    cost = float(parts[1])
    return (time, cost)

# # Variable parameters
# parameters = [
#     ("", ['-r']),
#     ('nnodes', [500]),
#     ('seed', list(range(0, 5))),
#     ('config', [0, 1])
# ]

# pp = generate_parameters(parameters)

# f = open(output_file, "w")

# for p in pp:
#     rawresult = run_experiment(p)
#     result = parse_result(rawresult)
#     f.write("{}, {}, {}\n".format(p, result[0], result[1]))

# f.close()

ninstances = 4
configs = [0, 1, 2]
f = open(output_file, "w")

f.write("{}".format(len(configs)))
for config in configs:
    f.write(",{}".format(config))
f.write("\n")

for n in range(ninstances):
    parameters = [
        ("", ['-r']),
        ('nnodes', [2000]),
        ('seed', [n*100]),
        ('config', configs)
    ]

    f.write("{}".format(n))
    pp = generate_parameters(parameters)
    for p in pp:
        rawresult = run_experiment(p)
        result = parse_result(rawresult)
        f.write(",{}".format(result[1]))
    f.write("\n")

f.close()
