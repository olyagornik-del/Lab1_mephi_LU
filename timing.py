import time


# запускаем функцию, возвращаем результат и время в секундах
def measure_time(func, *args):
    start = time.perf_counter()
    result = func(*args)
    elapsed = time.perf_counter() - start
    return result, elapsed
