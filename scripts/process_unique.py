

import sys


reps = []
prev_event = ""
content = {}
for line in sys.stdin:
    print(line.rstrip(),";unique_content=1.00;repeat_ranges=None(0-0)",sep="")

