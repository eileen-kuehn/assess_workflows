from __future__ import division
import math


def absolute_error(values):
    mean = sum(values)/len(values)
    try:
        std_error = math.sqrt(sum([(error-mean)**2 for error in values])/float(len(values)-1))
    except ZeroDivisionError:
        std_error = float("Inf")
    return mean, std_error


def relative_error(values):
    mean, std_error = absolute_error(values)
    return mean, std_error / mean


def errors(values):
    mean, std_error = absolute_error(values)
    return mean, std_error, std_error / mean


def uncorrelated_relative_error(values):
    result = 0
    mean = 0
    for left_value, right_value in values:
        mean += abs(left_value - right_value)
        try:
            result += abs(left_value - right_value) / (left_value + right_value)
        except ZeroDivisionError:
            pass
    return mean / len(values), result / len(values)
