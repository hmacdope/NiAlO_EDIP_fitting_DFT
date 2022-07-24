import numpy as np

with open("extracted_states.xyz", 'r') as f:
    lines = f.readlines()
    conflines = []
    break_next = False
    for line in reversed(lines):
        conflines.append(line)
        if break_next:
            break
        if "Energy:" in line:
            break_next = True
with open("final_state.xyz", 'w') as w:
    for line in reversed(conflines):
        w.write(line)
