""" Modules"""
import re
import os
import socket
import platform
import subprocess

def get_processor_name() -> str:
    """
        Get the name of the processor on the current system
    """
    if platform.system() == "Windows":
        return platform.processor()
    if platform.system() == "Darwin":
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + '/usr/sbin'
        cmd = subprocess.run(['sysctl','-n','machdep.cpu.brand_string'],
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=True)
        return cmd.stdout.strip()
    if platform.system() == "Linux":
        command = "grep model /proc/cpuinfo"
        all_info = subprocess.check_output(command, shell=True).decode('utf-8').strip()
        for line in all_info.split("\n"):
            if "model name" in line:
                return re.sub(".*model name.*:", "", line,1)
    return ""

def get_sys_info(logger = None) -> None:
    """
        Print information about the current system
    """
    try:
        if logger is not None:
            logger.info("=======================================================")
            logger.info("                   SYSTEM INFORMATION                  ")
            logger.info("=======================================================")
            logger.info('Platform:            %s (%s)', platform.system(),platform.release())
            logger.info('Architecture:        %s', platform.machine())
            logger.info('Hostname:            %s', socket.gethostname())
            logger.info('Processor:           %s', get_processor_name())
            logger.info("=======================================================")
    
        else:
            print("=======================================================")
            print("                   SYSTEM INFORMATION                  ")
            print("=======================================================")
            print('Platform:            %s (%s)', platform.system(),platform.release())
            print('Architecture:        %s', platform.machine())
            print('Hostname:            %s', socket.gethostname())
            print('Processor:           %s', get_processor_name())
            print("=======================================================")
    
    except:
        logger.info("=======================================================")
        logger.info("                   SYSTEM INFORMATION                  ")
        logger.info("=======================================================")
        logger.info('Platform:            %s (%s)', platform.system(),platform.release())
        logger.info('Architecture:        %s', platform.machine())
        logger.info('Hostname:            %s', socket.gethostname())
        logger.info('Processor:           %s   NAME NOT AVAILABLE')
        logger.info("=======================================================")
    
