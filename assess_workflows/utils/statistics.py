from __future__ import division
import math


def mean(values):
    return sum(values)/len(values)


def standard_deviation(values):
    """
    corrected standard variation for samples

    :param values:
    :return:
    """
    try:
        return math.sqrt(sum(_distances_to_mean(values)) / (len(values) - 1))
    except ZeroDivisionError:
        return values[0]


def population_standard_deviation(values):
    return math.sqrt(sum(_distances_to_mean(values)) / len(values))


def _distances_to_mean(values):
    mean = sum(values)/len(values)
    return [(value - mean)**2 for value in values]


def absolute_error(values):
    mean = sum(values)/len(values)
    try:
        # sample standard deviation
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
    """
    This approach follows the anomaly detection approach. Trying to normalise results looking from
    one side. But apparently, this is not correct.
    :param values:
    :return:
    """
    result = 0
    for left_distance, left_max_distance, right_distance, right_max_distance in values:
        result += (((abs(left_distance - right_distance) / 2) / left_max_distance) +
                   ((abs(left_distance - right_distance) / 2) / right_max_distance)) / 2
    return result / len(values)


def uncorrelated_relative_max_distance_deviation(values, expected_value=None):
    result = 0
    for left_distance, left_max_distance, right_distance, right_max_distance in values:
        if expected_value is not None:
            result += (abs(expected_value - left_distance) +
                           abs(expected_value - right_distance)) / 2 / \
                      (left_max_distance + right_max_distance)
        else:
            result += (abs(left_distance - right_distance) / 2) / \
                      (left_max_distance + right_max_distance)
    return result / len(values)


def uncorrelated_relative_deviation_and_standard_error(values, expected_value=None):
    """
    Method returns the mean relative distance derivation and its standard error.
    The standard error of the mean relative deviation estimates how far the sample mean is likely
    to be from the population mean. Thus the corrected sample standard deviation (n-1) is used for
    calculation.

    :param values:
    :return:
    """
    errors = []
    for left_distance, left_max_distance, right_distance, right_max_distance in values:
        if expected_value is not None:
            errors.append((abs(expected_value - left_distance) +
                           abs(expected_value - right_distance)) / 2 /
                          (left_max_distance + right_max_distance))
        else:
            errors.append((abs(left_distance - right_distance) / 2) /
                          (left_max_distance + right_max_distance))
    mean = sum(errors) / len(errors)
    sample_sd = math.sqrt(sum([(error - mean)**2 for error in errors]) / (len(errors) - 1)) / \
        math.sqrt(len(errors))
    return mean, sample_sd
