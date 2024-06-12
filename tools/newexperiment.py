#!/usr/bin/env python3

from subprocess import run, PIPE, Popen
import argparse
import json

parser = argparse.ArgumentParser(description="Run experiments")
parser.add_argument("--output", type=str, default="out.txt")
parser.add_argument("--executable", type=str, default="./main")
parser.add_argument("--timelimit", type=int, default=120)
parser.add_argument("--ninstances", type=int, default=5)
parser.add_argument("--nnodes", type=int, default=1000)
parser.add_argument("--parallel", type=int, default=1)
parser.add_argument("--costortime", type=str, default="cost", choices=["cost", "time"])
parser.add_argument("configs", type=str)

args = parser.parse_args()

print(args)

output_file = args.output
executable_name = args.executable
pfflag = "--parsefriendly"
timelimit = "{}".format(args.timelimit)
parallel_tasks = args.parallel
nnodes = args.nnodes
ninstances = args.ninstances
configs = json.loads(args.configs)
cost_or_time = 0 if args.costortime == "time" else 1

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

def parse_result(result: str):
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
        print("Running: {}".format(command))
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
