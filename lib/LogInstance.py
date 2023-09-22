"""
loguru 封装类，导入即可直接使用
# 当前文件名 logger.py
"""
import time
from functools import wraps
import os
import datetime
import loguru

# 单例类的装饰器
def singleton_class_decorator(cls):
    """
    装饰器，单例类的装饰器
    """
    # 在装饰器里定义一个字典，用来存放类的实例。
    _instance = {}

    # 装饰器，被装饰的类
    @wraps(cls)
    def wrapper_class(*args, **kwargs):
        # 判断，类实例不在类实例的字典里，就重新创建类实例
        if cls not in _instance:
            # 将新创建的类实例，存入到实例字典中
            _instance[cls] = cls(*args, **kwargs)
        # 如果实例字典中，存在类实例，直接取出返回类实例
        return _instance[cls]

    # 返回，装饰器中，被装饰的类函数
    return wrapper_class


@singleton_class_decorator
class Logger:

    def __init__(self, withDate, tid_target_filter="", _prep=""):
        self.withDate = withDate
        self.tid_target_filter = tid_target_filter
        self._prep = _prep
        self.logger_add()

    def get_project_path(self, project_path=None):
        if project_path is None:
            # 当前项目文件的，绝对真实路径
            # 路径，一个点代表当前目录，两个点代表当前目录的上级目录
            project_path = os.path.realpath((os.path.dirname(os.path.dirname(__file__))))

        # 返回当前项目路径
        return project_path

    def get_log_path(self):
        # 项目目录
        project_path = self.get_project_path()
        # 项目日志目录
        project_log_dir_temp = os.path.join(project_path, 'log')
        project_log_dir = os.path.join(project_log_dir_temp, self.withDate)
        # 日志文件名
        # project_log_filename = 'runtime_{}'.format(datetime.date.today())+"/"+self.tid_target_filter
        project_log_filename = self.tid_target_filter
        # 日志文件路径
        project_log_path = os.path.join(project_log_dir, project_log_filename)
        # 返回日志路径
        return project_log_path

    def logger_add(self):
        # self.logger.add(self.get_log_path(), level='DEBUG',
        #                 format='{time:YYYYMMDD HH:mm:ss} - '  # 时间
        #                        "{process.name} | "  # 进程名
        #                        "{thread.name} | "  # 进程名
        #                        '{module}.{function}:{line} - {level} -{message}',  # 模块名.方法名:行号
        #                 rotation="10 MB")

        # 错误日志 |{process.name} | {thread.name}
        loguru.logger.add(
            os.path.join(self.get_log_path(), "ERROR/{time:YYYY-MM-DD-HH}.log"),
            format="{time:YYYY-MM-DD at HH:mm:ss}  | {module}.{function}:{line}| {level} | {message}",
            filter=lambda x: True if x["level"].name == "ERROR" else False,
            rotation="00:00", retention=7, level='ERROR', encoding='utf-8'
        )
        # # 成功日志
        # loguru.logger.add(
        #     os.path.join(self.get_log_path(), "SUCCESS/{time:YYYY-MM-DD}"+self._prep+".log"),
        #     format="{time:YYYY-MM-DD at HH:mm:ss} |{process.name} | {thread.name} | {module}.{function}:{line}| {level} | {message}",
        #     filter=lambda x: True if x["level"].name == "SUCCESS" else False,
        #     rotation="00:00", retention=7, level='SUCCESS', encoding='utf-8',
        # )

        # Default日志
        loguru.logger.add(
            os.path.join(self.get_log_path(),
                         "Info/{time:YYYY-MM-DD-HH}_yprep.log"),
            format="{time:YYYY-MM-DD at HH:mm:ss}  | {module}.{function}:{line}|{level} | {message}",
            filter=lambda x: True if '[_yprep]' in x['message'] else False,
            rotation="00:00", retention=7, level='DEBUG', encoding='utf-8'
        )
        # Default日志``
        loguru.logger.add(
            os.path.join(self.get_log_path(),
                         "Info/{time:YYYY-MM-DD-HH:mm}_yrcep.log"),
            format="{time:YYYY-MM-DD at HH:mm:ss}  | {module}.{function}:{line}|{level} | {message}",
            filter=lambda x: True if '[_yrcep]' in x['message'] else False,
            rotation="00:00", retention=7, level='DEBUG', encoding='utf-8'
        )
        loguru.logger.add(
            os.path.join(self.get_log_path(),
                         "Info/{time:YYYY-MM-DD-HH}_yphost.log"),
            format="{time:YYYY-MM-DD at HH:mm:ss}  | {module}.{function}:{line}|{level} | {message}",
            filter=lambda x: True if '[_yphost]' in x['message'] else False,
            rotation="00:00", retention=7, level='DEBUG', encoding='utf-8'
        )
        loguru.logger.add(
            os.path.join(self.get_log_path(),
                         "Info/{time:YYYY-MM-DD-HH}_ycom.log"),
            format="{time:YYYY-MM-DD at HH:mm:ss}  | {module}.{function}:{line}|{level} | {message}",
            filter=lambda x: True if '[_ycom]' in x['message'] else False,
            rotation="00:00", retention=7, level='DEBUG', encoding='utf-8'
        )
        loguru.logger.add(
            os.path.join(self.get_log_path(),
                         "Info/{time:YYYY-MM-DD-HH}_common.log"),
            format="{time:YYYY-MM-DD at HH:mm:ss}  | {module}.{function}:{line}|{level} | {message}",
            filter=lambda x: False if '[_yphost]' in x['message'] or '[_yrcep]' in x['message'] or '[_yrcep]' in x['message'] or '[_ycom]' in x['message'] else True,
            rotation="00:00", retention=7, level='DEBUG', encoding='utf-8'
        )


    @property
    def get_logger(self):
        return loguru.logger

    @property
    def get_logger(self):
        return loguru.logger

'''
# 实例化日志类
'''



