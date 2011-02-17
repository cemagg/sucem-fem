# import hotshot
# from hotshot import stats

disc1f = DifferentialForm.OneformDiscretiser(mesh)
A = SystemMatrix.OneformOneform()(disc1f)

disc2f = DifferentialForm.TwoformDiscretiser(mesh)
# B = SystemMatrix.TwoformTwoform()(disc2f)

#AB = SystemMatrix.OneformTwoform()(disc1f, disc2f)

# def main():
#     A = SystemMatrix.OneformOneform()(disc1f)

# prof = hotshot.Profile("hotshot_oneformsysmat")
# prof.runcall(main)
# prof.close()


# s = stats.load("hotshot_oneformsysmat")
# s.sort_stats("time")
