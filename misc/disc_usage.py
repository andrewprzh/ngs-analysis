#!/usr/bin/env python3
import subprocess
import requests
import sys
import os
from os.path import abspath, dirname, realpath, join, isfile
import time
import datetime

# ---------- Configuration ----------
DISK_NAME = "abga/work"              # mount point or device name
THRESHOLD = 95               # % usage before sending alert
CHECK_INTERVAL = 4 * 3600    # 4 hours in seconds
WEEKLY_INTERVAL = 7 * 24 * 3600  # 1 week in seconds

MAILGUN_DOMAIN = open("/home/andreyp/.mg/.domain").readlines()[0].strip()   # e.g. sandboxXXX.mailgun.org
MAILGUN_API_KEY = open("/home/andreyp/.mg/.api").readlines()[0].strip()
FROM_EMAIL = "postmaster@sandbox04a9219346e6491a88f560e5d5db7551.mailgun.org"
TO_EMAILS = list(map(lambda x: x.strip(), open("/home/andreyp/.mg/.emails").readlines()))
# -----------------------------------

def get_disk_usage(disk_name):
    """Return usage percentage for given disk (mount or device name)."""
    try:
        result = subprocess.run(["df", "-h"], capture_output=True, text=True, check=True)
        for line in result.stdout.splitlines()[1:]:
            fields = line.split()
            if disk_name in fields[0] or disk_name in fields[-1]:
                usage_str = fields[4]  # e.g. "75%"
                return int(usage_str.strip("%")), line
        print(f"Disk {disk_name} not found in df output.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print("Error running df:", e)
        sys.exit(1)

def send_mailgun_email(subject, text):
    """Send an email using Mailgun API."""
    response = requests.post(
        f"https://api.mailgun.net/v3/{MAILGUN_DOMAIN}/messages",
        auth=("api", MAILGUN_API_KEY),
        data={
            "from": FROM_EMAIL,
            "to": TO_EMAILS,
            "subject": subject,
            "text": text,
        }
    )
    return response

def main():
    last_weekly_report = 0  # epoch timestamp of last weekly email

    while True:
        usage, df_line = get_disk_usage(DISK_NAME)
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{now}] Disk {DISK_NAME} usage: {usage}%")

        # --- Check threshold ---
        if usage >= THRESHOLD:
            subject = f"[ALERT] Disk usage on {DISK_NAME} is {usage}%"
            text = f"Warning: Disk {DISK_NAME} is at {usage}% capacity.\n\nDF line:\n{df_line}"
            resp = send_mailgun_email(subject, text)
            print("Alert sent:", resp.status_code, resp.text)

        # --- Weekly report ---
        now_ts = time.time()
        if now_ts - last_weekly_report >= WEEKLY_INTERVAL:
            subject = f"[REPORT] Weekly disk usage report for {DISK_NAME}"
            text = f"Weekly report:\nDisk {DISK_NAME} usage is {usage}%.\n\nDF line:\n{df_line}"
            resp = send_mailgun_email(subject, text)
            print("Weekly report sent:", resp.status_code, resp.text)
            last_weekly_report = now_ts

        # Wait until next check
        time.sleep(CHECK_INTERVAL)

if __name__ == "__main__":
    main()
