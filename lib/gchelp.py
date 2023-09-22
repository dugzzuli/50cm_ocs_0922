import datetime
from functools import wraps
import inspect
import time

import loguru
import sqlite3

def elapsed_time(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        file_name = inspect.getfile(func)
        loguru.logger.info(f"[ycom_]{func.__name__}, {file_name} 运行时间：{elapsed_time:.3f} 秒")
        
        return result
    
    return wrapper

def run_times_outt(func):
    count = 0
    def run_times_inn(*args,**kwargs):
        nonlocal count
        count += 1
        func(*args,**kwargs)
        print('函数func运行次数为：{}'.format(count))
    return run_times_inn
class funMy():
    @elapsed_time
    def my_func(self):
        time.sleep(1) # 模拟代码运行
@elapsed_time
def my_func():
    time.sleep(1) # 模拟代码运行


def __clear_env():
    
    for key in globals().keys():

        if not key.startswith("__"): # 排除系统内建函数

            globals().pop(key)

