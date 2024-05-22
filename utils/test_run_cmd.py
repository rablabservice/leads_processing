#!/usr/bin/env python

import argparse
import time


def test_run_cmd(max_duration=10, wait_time=2):
    start_time = time.time()
    while True:
        elapsed = time.time() - start_time
        if elapsed > max_duration:
            break
        else:
            print(f"{elapsed:.2f} seconds elapsed...")
            time.sleep(wait_time)
    print(f"Total time: {elapsed:.2f} seconds")


def _parse_args():
    parser = argparse.ArgumentParser(description="Test run command")
    parser.add_argument(
        "-m",
        "--max_duration",
        type=float,
        default=10,
        help="Maximum duration to run the test",
    )
    parser.add_argument(
        "-w",
        "--wait_time",
        type=float,
        default=2,
        help="Time to wait between each iteration",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    test_run_cmd(args.max_duration, args.wait_time)
