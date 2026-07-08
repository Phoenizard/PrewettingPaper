"""Timestamped, flushed progress logging for long-running loops."""

import time


def log(msg):
    print(time.strftime("[%H:%M:%S] ") + msg, flush=True)
