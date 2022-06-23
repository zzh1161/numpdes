###
 # @Author: zhang-zh 1563833973@qq.com
 # @Date: 2022-05-25 14:22:41
 # @LastEditors: zhang-zh 1563833973@qq.com
 # @LastEditTime: 2022-05-25 14:22:41
 # @FilePath: /undefined/home/zhang-zh/workspace/MyGithub/numPDEs/.zshrc
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 
# Proxy for Windows
export hostip=$(cat /etc/resolv.conf |grep -oP '(?<=nameserver\ ).*')
alias setss='export https_proxy="http://${hostip}:10808";export http_proxy="http://${hostip}:10808";export all_proxy="socks5://${hostip}:10808";'