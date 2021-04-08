from datetime import datetime


def closest_previous(arr, val):
    previous = []
    for i in range(len(arr)):
        if val >= arr[i]:
            previous.append(arr[i])
    return min(previous, key=lambda x: abs(x - val))


def log(str, indent=0):
    out = datetime.now().strftime("%H:%M:%S.%f") + (" " * 3 * (indent + 1)) + str
    print(out)
    with open("log_hydrosim.txt", "a") as file:
        file.write(out + "\n")


def error(str):
    out = datetime.now().strftime("%H:%M:%S.%f") + "   ERROR: " + str
    with open("log_hydrosim.txt", "a") as file:
        file.write(out + "\n")
    raise ValueError(str)
