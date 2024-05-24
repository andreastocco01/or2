#!/usr/bin/env python3

from subprocess import run, PIPE, Popen

output_file = "out.txt"
executable_name = "./main"
pfflag = "--parsefriendly"
timelimit = "60"
parallel_tasks = 1
nnodes = 1000
ninstances = 5
configs = [2, 12]
cost_or_time = 1 # 0 for time, 1 for cost 

def divide_list_into_blocks(lst, n):
    return [lst[i:i+n] for i in range(0, len(lst), n)]

def make_seed(n) -> int:
    return n*100

def build_command(parameters):
    seed = str(make_seed(parameters[0]))
    config = str(parameters[1])
    n = str(nnodes)
    command = [ executable_name,
                pfflag,
                "-r",
                "-t", timelimit,
                "--config", config,
                "--seed", seed,
                "--nnodes", n]
    return command

def parse_result(result: str) -> (float, float):
    result = result.replace("\n ", "")
    parts = result.split(";")
    time = float(parts[0])
    cost = float(parts[1])
    return (time, cost)

tasks = []

for n in range(ninstances):
    for config in configs:
        tasks.append((n, config,))

blocks = divide_list_into_blocks(tasks, parallel_tasks)

done = dict()

blockscount = len(blocks)
counter = 1

for block in blocks:
    print("Running block {}/{}".format(counter, blockscount));
    commands = []
    for task in block:
        command = build_command(task)
        commands.append(command)
    procs = [Popen(i, stdout=PIPE) for i in commands ]
    stdouts = [process.stdout for process in procs]
    for p in procs:
        p.wait()
    results = []
    for stdout in stdouts:
        last = next(stdout).decode('utf-8')
        parsed = parse_result(last)
        results.append(parsed)

    zipped = zip(block, results)
    for el in zipped:
        done[el[0]] = el[1]

    counter += 1


f = open(output_file, "w")

f.write("{}".format(len(configs)))
for config in configs:
    f.write(",{}".format(config))
f.write("\n")

for n in range(ninstances):
    f.write("{}".format(n))
    for config in configs:
        result = done[(n, config)][cost_or_time]
        f.write(",{}".format(result))
    f.write("\n")

f.close()
