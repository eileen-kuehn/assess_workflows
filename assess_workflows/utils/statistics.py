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
    for left_value, right_value in values:
        try:
            result += abs(left_value - right_value) / (left_value + right_value)
        except ZeroDivisionError:
            pass
    return result / len(values)


def uncorrelated_relative_distance_deviation(values):
    result = 0
    for left_distance, left_max_distance, right_distance, right_max_distance in values:
        result += (((abs(left_distance - right_distance) / 2) / left_max_distance) +
                   ((abs(left_distance - right_distance) / 2) / right_max_distance)) / 2
    return result / len(values)


def uncorrelated_relative_max_distance_deviation(values):
    result = 0
    for left_distance, left_max_distance, right_distance, right_max_distance in values:
        result += (abs(left_distance - right_distance) / 2) / \
                  (left_max_distance + right_max_distance)
    return result / len(values)
